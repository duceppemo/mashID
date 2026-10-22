"""Which references of a database belong to different organisms yet cannot be told apart.

Whether two organisms are separable is a property of the database, not of any one sample: the same
pair of references is the same distance apart in every run. Computing it once, at build time, gives
the per-sample annotation for free, and does so correctly whether or not Mash's winner-take-all
strategy is in use -- under winner-take-all the runner-up's identity is the count of hashes *unique*
to that reference, not a similarity, so the gap between the top two hits cannot be read as a margin.

The same table is the most direct evidence of a mislabelled reference: a genome named for one species
sitting a hair away from the references of another is either an inseparable pair or a bad label.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from mashid import MashIDError
from mashid.mash import dist_pairs
from mashid.metadata import NA, DbEntry, accession_from_query_id

log = logging.getLogger(__name__)

# Mash distance below which two references of different organisms are reported as inseparable. The
# default matches the default --ambiguity-margin: both ask "could this call just as well be the other".
DEFAULT_MAX_DISTANCE = 0.005
CONFUSABLE_SUFFIX = ".confusable.tsv"
COLUMNS = ["Accession", "Organism", "Confusable_Organism", "Distance", "Nearest_Accession"]

# Species placeholders: an unnamed reference makes no claim to a species, so it cannot contradict one.
PLACEHOLDER_SPECIES = {"sp.", "sp", "spp.", "bacterium", "unknown"}
MAX_NAMES_IN_NOTE = 3

# This note states a property of the database that holds for every sample hitting that reference, so
# it says nothing about the sample and does not count towards --fail-on note.
NOTE_PREFIX = "Not separable in this database from: "


def confusable_path(database: Path) -> Path:
    """``/x/db.msh`` -> ``/x/db.confusable.tsv``."""
    database = Path(database)
    stem = database.name[:-4] if database.name.endswith(".msh") else database.name
    return database.with_name(stem + CONFUSABLE_SUFFIX)


def _binomial(organism: str) -> tuple[str, str] | None:
    """(genus, species epithet) of a Latin binomial, or None if the name is not one."""
    tokens = organism.replace("[", "").replace("]", "").split()
    if tokens[:1] == ["Candidatus"]:
        tokens = tokens[1:]
    if len(tokens) < 2 or tokens[1] in PLACEHOLDER_SPECIES:
        return None
    return tokens[0], tokens[1]


def distinguishable(a: str, b: str) -> bool:
    """True when ``a`` and ``b`` are claims about genuinely different organisms.

    False for identical names, for names without a usable species (an unnamed "Genus sp." reference
    contradicts nothing), and for the same species epithet under two genera, which in practice means
    a renamed genus rather than two organisms (Mycobacterium abscessus / Mycobacteroides abscessus).
    """
    a, b = a.strip(), b.strip()
    if not a or not b or a in (NA, "unknown") or b in (NA, "unknown") or a == b:
        return False
    bin_a, bin_b = _binomial(a), _binomial(b)
    if bin_a is None or bin_b is None:
        return False
    # The same epithet under two genera is a renamed genus, not a second organism.
    return not (bin_a[1] == bin_b[1] and bin_a[0] != bin_b[0])


@dataclass(frozen=True)
class Neighbour:
    organism: str
    distance: float
    accession: str


def compute(database: Path, metadata: dict[str, DbEntry], threads: int = 1,
            max_distance: float = DEFAULT_MAX_DISTANCE,
            exe: str = "mash") -> dict[str, list[Neighbour]]:
    """All-pairs ``mash dist`` of the database against itself, keeping inseparable cross-organism pairs.

    Returns ``{accession: [Neighbour, ...]}``, nearest first, holding for each reference the closest
    reference of each *other* organism within ``max_distance``.
    """
    database = Path(database)
    if not database.is_file():
        raise MashIDError(f"Mash database not found: {database}")
    if not metadata:
        raise MashIDError(f"No reference metadata for {database}; write the sidecar first with "
                          f"make_mashID_db --annotate")

    def organism_of(query_id: str) -> tuple[str, str]:
        accession = accession_from_query_id(query_id)
        entry = metadata.get(accession)
        return accession, (entry.organism if entry is not None else NA)

    best: dict[str, dict[str, Neighbour]] = {}
    unnamed: set[str] = set()
    pairs = seen = 0
    for id_a, id_b, distance in dist_pairs(database, database, threads, max_distance, exe):
        seen += 1
        acc_a, org_a = organism_of(id_a)
        acc_b, org_b = organism_of(id_b)
        for acc, org in ((acc_a, org_a), (acc_b, org_b)):
            if org == NA:
                unnamed.add(acc)
        if acc_a == acc_b or not distinguishable(org_a, org_b):
            continue
        pairs += 1
        for acc, other_acc, other_org in ((acc_a, acc_b, org_b), (acc_b, acc_a, org_a)):
            current = best.setdefault(acc, {}).get(other_org)
            if current is None or distance < current.distance:
                best[acc][other_org] = Neighbour(other_org, distance, other_acc)
    if unnamed:
        # Sketch names that carry no accession the sidecar knows: every pair involving them is
        # skipped, so say so rather than quietly reporting that nothing is confusable.
        log.warning("%d reference(s) in the compared pairs have no organism in the metadata sidecar "
                    "and were skipped (for example %s); refresh it with make_mashID_db --annotate",
                    len(unnamed), ", ".join(sorted(unnamed)[:3]))
    log.info("Compared %d reference pair(s) within distance %g; %d cross-organism pair(s) over %d "
             "reference(s)", seen, max_distance, pairs, len(best))
    return {acc: sorted(n.values(), key=lambda x: (x.distance, x.organism))
            for acc, n in sorted(best.items())}


def write(path: Path, neighbours: dict[str, list[Neighbour]], metadata: dict[str, DbEntry]) -> Path:
    """Write the table, nearest pairs first so that mislabelled references sort to the top."""
    path = Path(path)
    rows = []
    for accession, items in neighbours.items():
        entry = metadata.get(accession)
        organism = entry.organism if entry is not None else NA
        rows.extend({"Accession": accession, "Organism": organism,
                     "Confusable_Organism": n.organism, "Distance": f"{n.distance:.6f}",
                     "Nearest_Accession": n.accession} for n in items)
    rows.sort(key=lambda r: (float(r["Distance"]), r["Accession"], r["Confusable_Organism"]))
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    log.info("Wrote %d row(s) for %d reference(s) to %s", len(rows), len(neighbours), path)
    return path


def read(path: Path) -> dict[str, list[Neighbour]]:
    """Read the table back into ``{accession: [Neighbour, ...]}``, nearest first."""
    path = Path(path)
    if not path.is_file():
        raise MashIDError(f"Confusable-reference table not found: {path}")
    out: dict[str, list[Neighbour]] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"Accession", "Confusable_Organism", "Distance"} - set(reader.fieldnames or [])
        if missing:
            raise MashIDError(f"{path}: missing column(s) {', '.join(sorted(missing))}")
        for row in reader:
            accession = (row.get("Accession") or "").strip()
            organism = (row.get("Confusable_Organism") or "").strip()
            if not accession or not organism:
                continue
            try:
                distance = float(row.get("Distance") or "")
            except ValueError:
                continue
            out.setdefault(accession, []).append(
                Neighbour(organism, distance, (row.get("Nearest_Accession") or "").strip()))
    for items in out.values():
        items.sort(key=lambda x: (x.distance, x.organism))
    return out


def load(database: Path, override: Path | None) -> dict[str, list[Neighbour]] | None:
    """``--db-confusable``, else the sidecar next to the database, else None (check disabled)."""
    if override is not None:
        table = read(override)
        log.info("Loaded inseparable-organism pairs for %d reference(s) from %s", len(table), override)
        return table
    path = confusable_path(database)
    if not path.is_file():
        return None
    table = read(path)
    log.info("Loaded inseparable-organism pairs for %d reference(s) from %s", len(table), path)
    if path.stat().st_mtime < Path(database).stat().st_mtime:
        log.warning("%s is older than the database; if the database was rebuilt, refresh it with "
                    "make_mashID_db --confusable", path.name)
    return table


def note_for(accession: str, organism: str, table: dict[str, list[Neighbour]] | None,
             exclude: set[str] | None = None) -> str | None:
    """The summary note for a top hit, or None when the database separates it from everything else."""
    if not table:
        return None
    skip = {organism} | (exclude or set())
    names: list[str] = []
    for neighbour in table.get(accession, []):
        if neighbour.organism not in skip and neighbour.organism not in names:
            names.append(neighbour.organism)
    if not names:
        return None
    shown = names[:MAX_NAMES_IN_NOTE]
    rest = len(names) - len(shown)
    return NOTE_PREFIX + ", ".join(shown) + (f" and {rest} more" if rest else "")
