# PlasmidNet

Interactive dashboard for exploring the [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) plasmid database (72,360+ plasmids).

## Features

- **Overview**: Topology, source, kingdom distributions; temporal trends; GC content and length distributions
- **Taxonomy Explorer**: Top 20 host genera visualization
- **AMR Analysis**: Antimicrobial resistance gene frequencies grouped by drug class
- **Plasmid Lookup**: Search individual plasmids by NCBI accession and view full metadata

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run locally
python app.py
```

Open http://localhost:8050 in your browser.

## Deploy to Heroku

```bash
heroku create plasmidnet
git push heroku main
```

## Deploy to AWS Elastic Beanstalk

```bash
pip install awsebcli
eb init -p python-3.12 plasmidnet
eb create plasmidnet-env
eb open
```

## Deploy to Render (free tier)

1. Connect your GitHub repo at [render.com](https://render.com)
2. Create a new Web Service
3. Set build command: `pip install -r requirements.txt`
4. Set start command: `gunicorn app:server --bind 0.0.0.0:$PORT`

## Tech Stack

- [Dash](https://dash.plotly.com/) + [Plotly](https://plotly.com/) for interactive visualizations
- Gunicorn for production serving


