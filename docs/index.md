# PlasmidNet Documentation

**PlasmidNet** is an interactive dashboard for exploring **182,362 bacterial plasmids** from [PLSDB 2025](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) and [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/). It provides comprehensive analytics for antimicrobial resistance, plasmid mobility, geographic distribution, and sequence analysis.

**Live Dashboard**: [plasmidnet-7302570f4818.herokuapp.com](https://plasmidnet-7302570f4818.herokuapp.com/)

**Source Code**: [github.com/lcerdeira/plasmidnet](https://github.com/lcerdeira/plasmidnet)

## Quick Start

```bash
pip install -r requirements.txt
python build_db.py    # downloads PLSDB data from Figshare, builds SQLite DB
python app.py         # open http://localhost:8050
```

## Contents

```{toctree}
:maxdepth: 2

features
api
seq-analysis
compare
analytics
database
deployment
```
