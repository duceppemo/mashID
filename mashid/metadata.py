"""Database metadata: accession parsing and the optional ``<db>.metadata.tsv`` sidecar.

The sidecar maps each reference accession to an organism name (and NCBI TaxID when known) so that
identification does not depend on parsing fasta headers. It is written by ``make_mashID_db`` and read by
``mashID`` when found next to the ``.msh`` file.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

from mashid import MashIDError
from mashid.mash import SketchInfo
from mashid.samples import strip_seq_extension
from mashid.taxonomy import clean_comment, organism_from_comment

# NCBI assembly accessions: GCF_000195955.2 / GCA_000195955.2
ACCESSION_RE = re.compile(r"(GC[AF]_\d{9}\.\d+)")

METADATA_COLUMNS = ["Accession", "Organism", "TaxID", "Length", "Hashes", "Source_File", "Description"]
NA = "NA"


@dataclass
class DbEntry:
    accession: str
    organism: str
    taxid: str = NA
    length: int | None = None  # reference sequence length in bp (from mash info)
    hashes: int | None = None  # sketch size actually stored for this reference
    source_file: str = ""
    description: str = ""


def accession_from_query_id(query_id: str) -> str:
    """'/db/fna/GCF_001667995.1_ASM166799v1_genomic.fna.gz' -> 'GCF_001667995.1'; otherwise the file stem."""
    name = Path(query_id).name
    match = ACCESSION_RE.search(name)
    if match:
        return match.group(1)
    return strip_seq_extension(name) or name


def sidecar_path(database: Path) -> Path:
    """'/db/foo.msh' -> '/db/foo.metadata.tsv'."""
    database = Path(database)
    stem = database.name[:-4] if database.name.endswith(".msh") else database.name
    return database.with_name(stem + ".metadata.tsv")


def _normalise_header(fieldnames: list[str] | None) -> dict[str, str]:
    """Map lowercase, punctuation-free column names to the actual header names."""
    if not fieldnames:
        return {}
    return {re.sub(r"[^a-z0-9]", "", f.lower()): f for f in fieldnames}


def _open_table(path: Path):
    text = path.read_text(encoding="utf-8-sig")
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    return csv.DictReader(text.splitlines(), delimiter=delimiter)


def read_metadata(path: Path) -> dict[str, DbEntry]:
    """Read a sidecar (or any TSV/CSV with accession/organism[/taxid] columns) into a dict by accession."""
    path = Path(path)
    if not path.is_file():
        raise MashIDError(f"Metadata file not found: {path}")
    reader = _open_table(path)
    cols = _normalise_header(reader.fieldnames)
    acc_col = cols.get("accession") or cols.get("assemblyaccession") or cols.get("acc")
    org_col = cols.get("organism") or cols.get("organismname") or cols.get("species") or cols.get("name")
    if not acc_col or not org_col:
        raise MashIDError(
            f"{path}: expected columns 'Accession' and 'Organism' (optionally 'TaxID'); found {reader.fieldnames}"
        )
    tax_col = cols.get("taxid") or cols.get("taxonomyid") or cols.get("organismtaxid")
    len_col = cols.get("length")
    hash_col = cols.get("hashes")
    src_col = cols.get("sourcefile")
    desc_col = cols.get("description")

    def _int(value: str | None) -> int | None:
        try:
            return int(value) if value not in (None, "", NA) else None
        except ValueError:
            return None
    entries: dict[str, DbEntry] = {}
    for row in reader:
        acc = (row.get(acc_col) or "").strip()
        if not acc:
            continue
        entries[acc] = DbEntry(
            accession=acc,
            organism=(row.get(org_col) or "").strip() or NA,
            taxid=(row.get(tax_col) or "").strip() or NA if tax_col else NA,
            length=_int(row.get(len_col)) if len_col else None,
            hashes=_int(row.get(hash_col)) if hash_col else None,
            source_file=(row.get(src_col) or "").strip() if src_col else "",
            description=(row.get(desc_col) or "").strip() if desc_col else "",
        )
    return entries


def write_metadata(path: Path, entries: list[DbEntry]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(METADATA_COLUMNS)
        for e in entries:
            writer.writerow([e.accession, e.organism, e.taxid or NA,
                             NA if e.length is None else e.length, NA if e.hashes is None else e.hashes,
                             e.source_file, e.description])


def read_ncbi_assembly_report(path: Path) -> dict[str, tuple[str, str]]:
    """Parse NCBI ``datasets`` ``assembly_data_report.jsonl`` into {accession: (organism, taxid)}."""
    path = Path(path)
    if not path.is_file():
        raise MashIDError(f"Assembly report not found: {path}")
    result: dict[str, tuple[str, str]] = {}
    with open(path, encoding="utf-8") as fh:
        for n, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise MashIDError(f"{path}: line {n} is not valid JSON: {exc}") from exc
            acc = record.get("accession") or record.get("assemblyInfo", {}).get("assemblyAccession")
            organism = record.get("organism", {}) or {}
            name = normalise_organism_name(organism.get("organismName") or NA)
            taxid = str(organism.get("taxId") or NA)
            if acc:
                result[acc] = (name, taxid)
                # Also register the paired accession (GCF <-> GCA) so either naming resolves.
                paired = record.get("pairedAccession")
                if paired:
                    result.setdefault(paired, (name, taxid))
    return result


def normalise_organism_name(name: str) -> str:
    """Strip strain-level text from an NCBI organism name.

    "Mycobacterium tuberculosis variant bovis BCG str. Sweden" -> "Mycobacterium tuberculosis variant bovis"
    "Mycobacterium sp. JS623" is kept as is. Names that do not parse are returned unchanged.
    """
    if not name or name == NA:
        return NA
    parsed = organism_from_comment(name, fallback="")
    return parsed or name


def lookup(entries: dict[str, DbEntry] | None, accession: str) -> DbEntry | None:
    if not entries:
        return None
    entry = entries.get(accession)
    if entry is None and "." in accession:  # allow version-less matches
        base = accession.rsplit(".", 1)[0]
        entry = next((e for a, e in entries.items() if a.rsplit(".", 1)[0] == base), None)
    return entry


def entries_from_sketches(sketches: list[SketchInfo],
                          annotations: dict[str, tuple[str, str]] | None = None) -> tuple[list[DbEntry], int]:
    """Build sidecar entries for every sketch of a database.

    Organism/TaxID come from ``annotations`` ({accession: (organism, taxid)}) when the accession is
    listed, else the organism is parsed from the sketch comment. Returns (entries, n_unannotated).
    """
    annotations = annotations or {}
    entries: list[DbEntry] = []
    unannotated = 0
    for sk in sketches:
        accession = accession_from_query_id(sk.query_id)
        if accession in annotations:
            organism, taxid = annotations[accession]
        else:
            organism, taxid = organism_from_comment(sk.comment, fallback=accession), NA
            unannotated += 1
        entries.append(DbEntry(accession, organism, taxid or NA, sk.length or None, sk.hashes or None,
                               sk.query_id, clean_comment(sk.comment)))
    return entries, unannotated
