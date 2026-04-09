#!/usr/bin/env bash
# =============================================================================
# PlasmidNet Pipeline for r6i.2xlarge (8 vCPU, 64GB RAM, /data = 1TB free)
# Run: nohup bash ec2_run.sh > pipeline.log 2>&1 &
# =============================================================================
set -euo pipefail

NCPU=8
WORKDIR="/data/plasmidnet_pipeline"
ACCESSION_FILE="$WORKDIR/ncbi_accessions.txt"
BATCH_SIZE=200
FASTA_DIR="$WORKDIR/fasta_batches"
RESULTS_DIR="$WORKDIR/results"

echo "============================================="
echo "PlasmidNet MOBsuite + AMRFinderPlus Pipeline"
echo "CPUs: $NCPU | RAM: $(free -h | awk '/Mem/{print $2}')"
echo "Start: $(date)"
echo "============================================="

mkdir -p "$FASTA_DIR" "$RESULTS_DIR"

# =============================================================================
# Step 1: Install Miniconda + MOBsuite + AMRFinderPlus
# =============================================================================
echo ""
echo ">>> Step 1: Installing dependencies..."

if [ ! -d "$HOME/miniconda3" ]; then
    echo "  Installing Miniconda..."
    curl -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
fi
export PATH="$HOME/miniconda3/bin:$PATH"
eval "$(conda shell.bash hook)"

if ! conda env list 2>/dev/null | grep -q plasmidnet; then
    echo "  Creating conda env with MOBsuite + AMRFinderPlus..."
    conda create -y -n plasmidnet -c bioconda -c conda-forge \
        mob_suite amrfinderplus
fi
conda activate plasmidnet

# GNU parallel
if ! command -v parallel &>/dev/null; then
    sudo yum install -y parallel 2>/dev/null || conda install -y -c conda-forge parallel
fi

# AMRFinderPlus database
if [ ! -d "$WORKDIR/amrfinder_db" ]; then
    echo "  Downloading AMRFinderPlus database..."
    amrfinder_update -d "$WORKDIR/amrfinder_db"
fi

echo "  mob_typer: $(mob_typer --version 2>&1 | head -1)"
echo "  amrfinder: $(amrfinder --version 2>&1 | head -1)"

# =============================================================================
# Step 2: Download sequences from NCBI in batches
# =============================================================================
echo ""
echo ">>> Step 2: Downloading sequences..."

TOTAL_ACC=$(wc -l < "$ACCESSION_FILE")
echo "  Accessions: $TOTAL_ACC"

# Split into batches
mkdir -p "$WORKDIR/acc_batches"
split -l $BATCH_SIZE "$ACCESSION_FILE" "$WORKDIR/acc_batches/batch_"
NBATCHES=$(ls "$WORKDIR/acc_batches"/batch_* | wc -l)
echo "  Batches: $NBATCHES (${BATCH_SIZE} each)"

download_batch() {
    local batch_file="$1"
    local out_file="$2"
    [ -s "$out_file" ] && return 0  # skip if already downloaded
    local ids
    ids=$(tr '\n' ',' < "$batch_file")
    ids=${ids%,}
    for attempt in 1 2 3 4 5; do
        if curl -sf --max-time 120 \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ids}&rettype=fasta&retmode=text" \
            -o "$out_file" 2>/dev/null; then
            if [ -s "$out_file" ] && head -1 "$out_file" | grep -q "^>"; then
                return 0
            fi
        fi
        sleep $((attempt * 3))
    done
    echo "WARN: failed $(basename "$batch_file")" >&2
    return 0
}
export -f download_batch

DONE_DL=0
for bf in "$WORKDIR/acc_batches"/batch_*; do
    out="$FASTA_DIR/$(basename "$bf").fasta"
    if [ -s "$out" ]; then
        DONE_DL=$((DONE_DL + 1))
    fi
done
echo "  Already downloaded: $DONE_DL / $NBATCHES"

# Download 3 concurrent (NCBI rate limit)
ls "$WORKDIR/acc_batches"/batch_* | while read bf; do
    echo "$bf $FASTA_DIR/$(basename "$bf").fasta"
done | parallel -j 3 --colsep ' ' --bar download_batch {1} {2}

SEQ_COUNT=$(grep -c "^>" "$FASTA_DIR"/*.fasta 2>/dev/null | tail -1 | cut -d: -f2 || \
            cat "$FASTA_DIR"/*.fasta | grep -c "^>" || echo 0)
echo "  Downloaded: $SEQ_COUNT sequences"

# =============================================================================
# Step 3: Run MOBsuite mob_typer
# =============================================================================
echo ""
echo ">>> Step 3: Running MOBsuite ($NCPU parallel)..."

MOB_OUTDIR="$WORKDIR/mob_out"
mkdir -p "$MOB_OUTDIR"

# Process each batch FASTA with mob_typer
run_mob_batch() {
    local fasta="$1"
    local outfile="$2"
    [ -s "$outfile" ] && return 0  # skip if done
    mob_typer --infile "$fasta" --out_file "$outfile" 2>/dev/null || true
}
export -f run_mob_batch

ls "$FASTA_DIR"/*.fasta | while read f; do
    echo "$f $MOB_OUTDIR/$(basename "$f" .fasta)_mob.txt"
done | parallel -j "$NCPU" --colsep ' ' --bar run_mob_batch {1} {2}

# Combine results
echo "  Combining MOBsuite results..."
first_file=$(ls "$MOB_OUTDIR"/*_mob.txt 2>/dev/null | head -1)
if [ -n "$first_file" ] && [ -f "$first_file" ]; then
    head -1 "$first_file" > "$RESULTS_DIR/mob_results.csv"
    for f in "$MOB_OUTDIR"/*_mob.txt; do
        tail -n +2 "$f" >> "$RESULTS_DIR/mob_results.csv" 2>/dev/null
    done
fi
MOB_LINES=$(wc -l < "$RESULTS_DIR/mob_results.csv" 2>/dev/null || echo 1)
echo "  MOBsuite: $((MOB_LINES - 1)) plasmids typed"

# =============================================================================
# Step 4: Run AMRFinderPlus
# =============================================================================
echo ""
echo ">>> Step 4: Running AMRFinderPlus ($NCPU parallel)..."

AMR_OUTDIR="$WORKDIR/amr_out"
mkdir -p "$AMR_OUTDIR"

run_amr_batch() {
    local fasta="$1"
    local outfile="$2"
    [ -s "$outfile" ] && return 0
    amrfinder -n "$fasta" -d "$WORKDIR/amrfinder_db/latest" \
        -o "$outfile" --plus 2>/dev/null || true
}
export -f run_amr_batch
export WORKDIR

ls "$FASTA_DIR"/*.fasta | while read f; do
    echo "$f $AMR_OUTDIR/$(basename "$f" .fasta)_amr.tsv"
done | parallel -j "$NCPU" --colsep ' ' --bar run_amr_batch {1} {2}

# Combine results
echo "  Combining AMRFinderPlus results..."
first_file=$(ls "$AMR_OUTDIR"/*_amr.tsv 2>/dev/null | head -1)
if [ -n "$first_file" ] && [ -f "$first_file" ]; then
    head -1 "$first_file" > "$RESULTS_DIR/amr_results.tsv"
    for f in "$AMR_OUTDIR"/*_amr.tsv; do
        tail -n +2 "$f" >> "$RESULTS_DIR/amr_results.tsv" 2>/dev/null
    done
fi
AMR_LINES=$(wc -l < "$RESULTS_DIR/amr_results.tsv" 2>/dev/null || echo 1)
echo "  AMRFinderPlus: $((AMR_LINES - 1)) annotations"

# =============================================================================
# Step 5: Summary
# =============================================================================
echo ""
echo "============================================="
echo "PIPELINE COMPLETE"
echo "End: $(date)"
echo "============================================="
echo ""
echo "Results in: $RESULTS_DIR/"
ls -lh "$RESULTS_DIR/"
echo ""
echo "Download from your Mac:"
echo "  scp -i ~/.ssh/dragon.pem ec2-user@44.211.173.63:$RESULTS_DIR/mob_results.csv ~/GitHub/PlasmidNet/"
echo "  scp -i ~/.ssh/dragon.pem ec2-user@44.211.173.63:$RESULTS_DIR/amr_results.tsv ~/GitHub/PlasmidNet/"
echo ""
echo "Then run:"
echo "  cd ~/GitHub/PlasmidNet"
echo "  python import_mobsuite_results.py mob_results.csv amr_results.tsv"
