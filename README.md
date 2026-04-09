[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19478464.svg)](https://doi.org/10.5281/zenodo.19478464)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-blue.svg)](https://www.rust-lang.org/)
![PlasmidNet](assets/logo.png)

Interactive dashboard for exploring the [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) plasmid database — all **72,556 plasmids** powered by a local SQLite database. No API rate limits, instant startup.





## Features

### Overview
Topology, source, kingdom distributions, temporal trends, GC content, and plasmid length distributions across the entire database.

### Taxonomy Explorer
Top 20 host genera with interactive search by genus or species — results from the full database in milliseconds.

### AMR Analysis
163 antimicrobial resistance drug classes and individual gene frequencies, grouped by drug class (beta-lactams, aminoglycosides, sulfonamides, etc.).

### Plasmid Viewer
Interactive circular plasmid map with zoom and hover tooltips. Enter any NCBI accession to visualize:
- **CDS genes** (forward/reverse strand) with product descriptions
- **AMR genes** highlighted in red (blaKPC-2, blaTEM, etc.)
- **MOB elements** (mobilization genes: mobA, traI, etc.)
- **Toxin-Antitoxin systems** (relE, higA, ccdA/B, vapB, etc.)
- **Virulence factors** (iutA, sitA, hlyA, etc.)
- **QAC resistance** (qacE, smr, etc.)
- **GC content** deviation ring (teal = above mean, red = below)

Works for both PLSDB plasmids and any NCBI GenBank accession (automatic fallback).

### Inc Groups & Mobility
Incompatibility group (replicon type) distribution, predicted mobility (conjugative/mobilizable/non-mobilizable), relaxase types, and MPF types.

### Correlations (n = 72,556)
Cross-feature co-occurrence analysis using the **entire database**:
1. **AMR Drug Class vs Inc Groups & Mobility** — heatmap and grouped bar
2. **Heavy Metal Resistance** — mercury, arsenic, copper, silver, tellurium gene frequencies and Inc group associations
3. **Phage & Mobile Elements** — distribution and mobility correlation
4. **pMLST Distribution** — scheme frequency, alleles, mobility
5. **Virulence Factors** (6,145 plasmids) — top genes, Inc group and mobility associations
6. **Toxin-Antitoxin Systems** — TA gene frequency and plasmid mobility
7. **QAC Resistance** (7,023 plasmids) — quaternary ammonium compound resistance vs Inc groups

### Plasmid Lookup
Search any accession to view metadata: length, GC, topology, taxonomy, source, creation date.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Build the SQLite database (one-time, ~30 seconds)
python build_db.py

# Run locally
python app.py
```

Open <http://localhost:8050> in your browser.

## Database

The dashboard uses a local SQLite database (`data/plsdb.db`, ~108 MB) built from [PLSDB Figshare](https://doi.org/10.6084/m9.figshare.27252609.v2) CSV exports. This eliminates all API calls and enables full-database analysis.

To rebuild after a PLSDB update:
```bash
# Download fresh CSVs to data/raw/ (see build_db.py)
python build_db.py
```

## Deploy to Heroku

```bash
heroku create plasmidnet
git config http.postBuffer 524288000  # for the 108 MB DB file
git push heroku main
heroku open
```

Note: Heroku has an ephemeral filesystem — the SQLite DB ships with git and is read-only at runtime, which works perfectly.

## Deploy to Render (free tier)

1. Connect your GitHub repo at [render.com](https://render.com)
2. Create a new Web Service
3. Set build command: `pip install -r requirements.txt`
4. Set start command: `gunicorn app:server --bind 0.0.0.0:$PORT`

## Tech Stack

- [Dash](https://dash.plotly.com/) + [Plotly](https://plotly.com/) for interactive visualizations
- SQLite for local full-database queries (72,556 plasmids, 251K AMR annotations)
- [PLSDB 2025](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) as data source
- NCBI GenBank as fallback for plasmid viewer (CDS + GC content)
- Gunicorn for production serving

<!-- ## Data Source

All data from PLSDB 2025 (Plasmid Database), version 2024_05_31_v2, maintained by the Computational Biology group at Saarland University. Plasmid viewer fetches CDS annotations and sequences from NCBI GenBank on demand. -->
