# PlasmidNet: Interactive Dashboard for 182,362 Bacterial Plasmids

## Description

PlasmidNet is an open-source interactive dashboard for exploring bacterial plasmid diversity, antimicrobial resistance, mobility, and geographic distribution. It integrates data from PLSDB 2025 (72,556 plasmids) and NCBI Nucleotide (109,806 additional plasmids) into a unified SQLite database with comprehensive analytics.

**Live Dashboard**: https://plasmidnet-7302570f4818.herokuapp.com/

**Source Code**: https://github.com/lcerdeira/plasmidnet

**Documentation**: https://plasmidnet.readthedocs.io/

## Files in this dataset

### Database

| File | Size | Description |
|------|------|-------------|
| `plsdb.db` | ~280 MB | SQLite database with all tables (nuccore, taxonomy, typing, amr, plasmidfinder, typing_markers, pgap_features, biosample_location, ncbi_extra) |

### Data Tables (CSV/TSV exports)

| File | Rows | Description |
|------|------|-------------|
| `plasmidnet_all_plasmids.csv` | 72,557 | All PLSDB plasmids with metadata (accession, length, GC, topology, taxonomy) |
| `plasmidnet_typing.csv` | 72,557 | MOB typing: Inc groups, mobility, relaxase, MPF, pMLST |
| `plasmidnet_amr.tsv` | 600,897 | AMR annotations from PLSDB + NCBI (gene, drug class, positions) |
| `plasmidnet_ncbi_extra.csv` | 109,807 | NCBI plasmids not in PLSDB (accession, length, organism, dates) |
| `plasmidnet_geography.csv` | 57,694 | Geographic data (country, lat/lng, host category, mobility) |

### Analysis Results

| File | Description |
|------|-------------|
| `correlations_cache.json` | Pre-computed correlation matrices (AMR vs Inc, heavy metals, VIR, TA, QAC) |
| `umap_embedding.npz` | UMAP 2D embedding of 15,000 plasmids (NumPy compressed) |
| `baseline_stats.json` | PLSDB baseline statistics for sequence analysis scoring |
| `mob_results.csv` | MOBsuite batch typing results for 109K NCBI plasmids |
| `amr_results.tsv` | AMRFinderPlus results for 109K NCBI plasmids (349,758 annotations) |

## Data Sources

- **PLSDB 2025** (version 2024_05_31_v2): https://doi.org/10.6084/m9.figshare.27252609.v2
- **NCBI Nucleotide**: Complete plasmid sequences harvested via E-utilities
- **PGAP annotations**: Filtered from PLSDB proteins.csv (TA, phage, transposase: 613,902 entries)
- **BioSample locations**: Geographic coordinates and host metadata from PLSDB biosample.csv

## Database Schema

```
nuccore (72,556 rows)
  NUCCORE_ACC, NUCCORE_Description, NUCCORE_Length, NUCCORE_GC,
  NUCCORE_Topology, NUCCORE_CreateDate, NUCCORE_Source,
  BIOSAMPLE_UID, TAXONOMY_UID

taxonomy (7,439 rows)
  TAXONOMY_UID, TAXONOMY_superkingdom, TAXONOMY_phylum, TAXONOMY_class,
  TAXONOMY_order, TAXONOMY_family, TAXONOMY_genus, TAXONOMY_species

typing (72,556 rows)
  NUCCORE_ACC, rep_type, relaxase_type, mpf_type, orit_type,
  predicted_mobility, predicted_host_range_overall_name,
  PMLST_scheme, PMLST_alleles, associated_pmid

amr (600,896 rows)
  NUCCORE_ACC, gene_symbol, gene_name, drug_class, antimicrobial_agent,
  input_gene_start, input_gene_stop, strand_orientation,
  sequence_identity, coverage_percentage

pgap_features (613,902 rows)
  NUCCORE_ACC, gene, product, locus_tag, category (TA/PHAGE/TRANSPOSASE),
  start, end, strand

biosample_location (23,522 rows)
  BIOSAMPLE_UID, lat, lng, country, ecosystem, host_category

ncbi_extra (109,806 rows)
  accession, description, length, topology, create_date, organism, taxid
```

## How to Use

### Quick start

```bash
# Clone and run
git clone https://github.com/lcerdeira/plasmidnet.git
cd plasmidnet
pip install -r requirements.txt
# Place plsdb.db in data/ directory
python app.py
# Open http://localhost:8050
```

### REST API

```python
import requests
base = "https://plasmidnet-7302570f4818.herokuapp.com/api"

# Search for KPC plasmids
kpc = requests.get(f"{base}/amr", params={"gene": "blaKPC"}).json()

# Get plasmid details
p = requests.get(f"{base}/plasmid/NZ_CP031107.1").json()
```

## Citation

If you use PlasmidNet in your research, please cite:

Cerdeira, L. (2026). PlasmidNet: Interactive Dashboard for Bacterial Plasmid Analytics. https://github.com/lcerdeira/plasmidnet

## License

MIT License
