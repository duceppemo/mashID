# Comparison with other tools

mashID answers one question quickly: which organism, among those in a database, is this isolate or
assembly? Other tools answer neighbouring questions. Use the one that matches yours.

| Tool | Method | Best at | Input | Resources | Where mashID differs |
| --- | --- | --- | --- | --- | --- |
| **mashID** | Mash screen (MinHash containment) against a sketch database | Fast species-level identification of isolates, with mixture and coverage flags | Reads or assemblies | Seconds per sample, < 2 GB | Curated names/TaxIDs per database, notes, run provenance, pipeline outputs. |
| [Mash screen](https://github.com/marbl/Mash) | The same containment estimate | The primitive that mashID wraps | Reads or assemblies | Same | mashID adds sample handling, database management, naming, notes and reports. |
| [sourmash](https://sourmash.readthedocs.io/) | FracMinHash sketches; `gather` decomposes a sample into references | Metagenome composition, large public databases (GTDB, GenBank) | Reads or assemblies | Minutes; large databases | sourmash `gather` reports every genome in a mixture with abundances; mashID reports the best hit and flags mixtures. |
| [Kraken2](https://github.com/DerrickWood/kraken2) / Bracken | Exact k-mer to LCA classification of each read | Read-level taxonomic profiles of metagenomes | Reads | Database in RAM (tens of GB for standard) | Per-read classification and full profiles; heavier and slower to set up. mashID has no per-read output. |
| [refseq_masher](https://github.com/phac-nml/refseq_masher) | Mash against a RefSeq sketch | Same idea as mashID against RefSeq | Reads or assemblies | Light | Fixed database; mashID lets you build and check your own and adds notes and provenance. |
| [KmerFinder](https://cge.food.dtu.dk/services/KmerFinder/) | k-mer overlap against a curated database | Web-based identification | Reads or assemblies | Web / CGE tools | Web service; mashID runs locally and offline. |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | Marker-gene placement in the GTDB reference tree, plus ANI | Authoritative taxonomy of assemblies, novel species | Assemblies only | Hours, > 60 GB RAM, 100 GB database | The most rigorous classification; not for reads or quick checks. |
| [FastANI](https://github.com/ParBLiSS/FastANI) / skani | Average nucleotide identity between genomes | Confirming species boundaries between two assemblies | Assemblies | Light | Pairwise ANI, not a database search; a good second step after mashID. |
| MLST / cgMLST, SNP typing | Allele or SNP calling against a scheme or reference | Within-species typing, lineage, outbreak analysis | Reads or assemblies | Varies | Resolution below the species level, where mashID stops. |

## When to use what

- **Isolate WGS, "what did we sequence?"** mashID, with a database for the expected group or a broad
  one. Confirm surprising results with FastANI against the top reference's assembly.
- **Metagenome or complex mixture:** sourmash `gather` or Kraken2 with Bracken; mashID only tells you
  that a mixture exists and names the strongest components.
- **Publication-grade taxonomy of a new assembly:** GTDB-Tk.
- **Lineage, variant, outbreak:** MLST, cgMLST or SNP-based tools after mashID has confirmed the
  species.
