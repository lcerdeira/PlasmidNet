# Database

PlasmidNet uses a local SQLite database built from PLSDB and NCBI data.

## Data Sources

| Source | Plasmids | Description |
|--------|----------|-------------|
| [PLSDB 2025](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) | 72,556 | Curated plasmid database (version 2024_05_31_v2) |
| [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) | 109,806 | Additional complete plasmids not in PLSDB |
| **Total** | **182,362** | |

## Tables

| Table | Rows | Description |
|-------|------|-------------|
| `nuccore` | 72,556 | Plasmid metadata (accession, length, GC, topology, dates) |
| `taxonomy` | 7,439 | Host taxonomy (kingdom to strain) |
| `typing` | 72,556 | MOB typing (Inc groups, mobility, relaxase, MPF, pMLST) |
| `amr` | 600,896 | AMR gene annotations from PLSDB + NCBI (gene, drug class, positions) |
| `plasmidfinder` | 60,619 | PlasmidFinder typing results |
| `typing_markers` | 333,648 | MOB suite detailed markers |
| `pgap_features` | 613,902 | PGAP annotations (TA, phage, transposase) |
| `biosample_location` | 23,522 | Geographic locations with host category |
| `ncbi_extra` | 109,806 | NCBI plasmids not in PLSDB |

## Building the Database

```bash
# 1. Download PLSDB CSVs from Figshare
mkdir -p data/raw
# Files are downloaded automatically by build_db.py

# 2. Build SQLite database
python build_db.py

# 3. (Optional) Harvest additional plasmids from NCBI
python ncbi_harvest.py

# 4. Rebuild with NCBI data
python build_db.py
```

## Figshare Files Used

Downloaded from [doi:10.6084/m9.figshare.27252609.v2](https://doi.org/10.6084/m9.figshare.27252609.v2):

| File | Size | Content |
|------|------|---------|
| nuccore.csv | 16 MB | Plasmid metadata |
| typing.csv | 22 MB | MOB typing results |
| amr.tsv | 36 MB | AMR annotations |
| taxonomy.csv | 2.3 MB | Host taxonomy |
| biosample.csv | 11 MB | Sample metadata + locations |
| plasmidfinder.csv | 15 MB | PlasmidFinder results |
| typing_markers.csv | 75 MB | MOB detail markers |
| proteins.csv | 637 MB | PGAP annotations (filtered to 60 MB for TA/phage/transposase) |

## Updating

To update when PLSDB releases a new version:

1. Download new CSVs from Figshare
2. Place in `data/raw/`
3. Run `python build_db.py`
4. Run `python ncbi_harvest.py` to catch new NCBI plasmids
5. Delete `data/correlations_cache.json` to rebuild correlations

## Database Size

- **SQLite file**: ~200 MB
- **Startup memory**: ~107 MB (correlations cached to JSON)
- **Heroku compatible**: Built during deploy via `bin/post_compile`
