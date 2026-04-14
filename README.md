[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19478464.svg)](https://doi.org/10.5281/zenodo.19478464)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-blue.svg)](https://www.python.org/)
![PlasmidNet](assets/logo.png)

Interactive dashboard for exploring **182,362 complete plasmids** from [PLSDB 2025](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) (72,556) + [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) (109,806 additional), powered by a local SQLite database. No API rate limits, instant startup.

[Documentation](https://plasmidnet.readthedocs.io/en/latest/database.html)
## Features

### Overview
Topology, source, kingdom distributions, temporal trends, GC content, and plasmid length distributions across the entire database.

### Taxonomy Explorer
Top 20 host genera with interactive search by genus or species.

### AMR Analysis
163 antimicrobial resistance drug classes and individual gene frequencies coloured by their real drug class from the database.

### Plasmid Viewer
Interactive circular plasmid map with zoom and hover tooltips. Enter any NCBI accession to visualize:
- **CDS genes** (forward/reverse strand) with product descriptions
- **AMR genes** highlighted in red (blaKPC-2, blaTEM, etc.)
- **MOB elements** (mobilization genes: mobA, traI, etc.)
- **Toxin-Antitoxin systems** (relE, higA, ccdA/B, vapB, etc.)
- **Virulence factors** (iutA, sitA, hlyA, etc.)
- **QAC resistance** (qacE, smr, etc.)
- **GC content** deviation ring

Works for both PLSDB plasmids and any NCBI GenBank accession (automatic fallback).

### Inc Groups & Mobility
Incompatibility group (replicon type) distribution, predicted mobility, relaxase types, and MPF types.

### Correlations (n = 72,556)
Cross-feature co-occurrence analysis using the **entire database**:
1. **AMR Drug Class vs Inc Groups & Mobility** -- heatmap and grouped bar
2. **Heavy Metal Resistance** -- mercury, arsenic, copper, silver, tellurium genes and Inc group associations
3. **Phage & Mobile Elements** -- distribution and mobility correlation
4. **pMLST Distribution** -- scheme frequency, alleles, mobility
5. **Virulence Factors** (6,145 plasmids) -- top genes, Inc group and mobility
6. **Toxin-Antitoxin Systems** (40,759 plasmids) -- vapB, ccdA/B, pemI/K from PGAP annotations
7. **QAC Resistance** (7,023 plasmids) -- quaternary ammonium compound resistance vs Inc groups

### Geography (57,751 geolocated plasmids)
1. **Global Plasmid Map** -- scatter map coloured by mobility, hover for accession/country/Inc group
2. **Country Comparison** -- mobility and Inc group distribution by top 20 countries
3. **Temporal Spread Animations** -- animated timeline (2000-2024) showing cumulative plasmid reports per country; separate animations for specific Inc groups (IncFIB, IncFII, etc.) and mobility types (conjugative, mobilizable, non-mobilizable)
4. **Plasmid Host & Source Correlations** -- plasmid distribution by host category (Human 11,464 / Animal 5,925 / Soil 1,576 / Water 1,332 / Food 1,059 / Environment 959), with heatmaps for Inc groups and AMR drug classes by host

### Plasmid Lookup
Search any accession to view metadata: length, GC, topology, taxonomy, source, creation date.

## Quick Start

```bash
pip install -r requirements.txt
python build_db.py    # downloads PLSDB CSVs from Figshare, builds SQLite DB
python app.py         # open http://localhost:8050
```

## Database

The dashboard uses a local SQLite database (`data/plsdb.db`, ~200 MB) built from [PLSDB Figshare](https://doi.org/10.6084/m9.figshare.27252609.v2) CSV exports including PGAP protein annotations. This eliminates all API calls and enables full-database analysis.

Tables: nuccore (72,556), taxonomy (7,439), typing (72,556), amr (251,138), plasmidfinder (60,619), typing_markers (333,648), pgap_features (613,902), biosample_location (23,522).

## Deploy to Heroku

```bash
heroku create plasmidnet
git push heroku main
heroku open
```

The `bin/post_compile` hook automatically downloads PLSDB CSVs from Figshare and builds the SQLite DB during deploy. No large files in git.

## Deploy to Render (free tier)

1. Connect your GitHub repo at [render.com](https://render.com)
2. Create a new Web Service
3. Set build command: `pip install -r requirements.txt && bash bin/post_compile`
4. Set start command: `gunicorn app:server --bind 0.0.0.0:$PORT`

## Tech Stack

- [Dash](https://dash.plotly.com/) + [Plotly](https://plotly.com/) for interactive visualizations
- SQLite for local full-database queries (72,556 plasmids, 251K AMR, 614K PGAP annotations)
- [PLSDB 2025](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) as data source (version 2024_05_31_v2)
- NCBI GenBank as fallback for plasmid viewer (CDS + GC content)
- Gunicorn for production serving
