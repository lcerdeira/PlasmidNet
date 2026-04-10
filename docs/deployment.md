# Deployment

## Local Development

```bash
git clone https://github.com/lcerdeira/plasmidnet.git
cd plasmidnet
pip install -r requirements.txt
python build_db.py    # ~2 min, downloads 200 MB from Figshare
python app.py         # http://localhost:8050
```

### Requirements

- Python 3.9+
- BLAST+ (for plasmid comparison): `brew install blast` or `apt install ncbi-blast+`
- ~500 MB disk space for the database

## Heroku

PlasmidNet is deployed on Heroku with automatic database building.

```bash
heroku create plasmidnet
git push heroku main
heroku open
```

### How it works

1. **heroku-community/apt** buildpack installs `ncbi-blast+` from `Aptfile`
2. **heroku/python** buildpack installs pip packages from `requirements.txt`
3. **bin/post_compile** hook:
   - Downloads PLSDB CSVs from Figshare
   - Filters proteins.csv for TA/phage/transposase genes
   - Harvests additional plasmids from NCBI
   - Builds the SQLite database
   - Pre-builds the correlations cache

### Heroku Limits

- **RAM**: 512 MB (Eco dyno). App uses ~107 MB at startup.
- **Slug size**: 500 MB limit. Current: ~140 MB.
- **Ephemeral filesystem**: Database is rebuilt on each deploy.

### Buildpacks

```
1. heroku-community/apt      (installs BLAST+)
2. heroku/python              (installs pip packages)
```

## Render (free tier)

1. Connect your GitHub repo at [render.com](https://render.com)
2. Create a new Web Service
3. Set build command: `pip install -r requirements.txt && bash bin/post_compile`
4. Set start command: `gunicorn app:server --bind 0.0.0.0:$PORT`

## Docker (optional)

```dockerfile
FROM python:3.12-slim
RUN apt-get update && apt-get install -y ncbi-blast+ curl && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
RUN bash bin/post_compile
CMD gunicorn app:server --bind 0.0.0.0:$PORT
```

## EC2 Pipeline (MOBsuite + AMRFinderPlus)

For typing the 109K NCBI plasmids:

```bash
# On EC2 (c6i.4xlarge or larger, 50 GB disk)
scp ec2_run.sh ec2-user@<ip>:/data/plasmidnet_pipeline/
scp data/ncbi_accessions.txt ec2-user@<ip>:/data/plasmidnet_pipeline/
ssh ec2-user@<ip>
cd /data/plasmidnet_pipeline
nohup bash run_now.sh > pipeline.log 2>&1 &

# Monitor
tail -f pipeline.log

# When done, download results
scp ec2-user@<ip>:/data/plasmidnet_pipeline/results/mob_results.csv .
scp ec2-user@<ip>:/data/plasmidnet_pipeline/results/amr_results.tsv .

# Import
python import_mobsuite_results.py mob_results.csv amr_results.tsv
```
