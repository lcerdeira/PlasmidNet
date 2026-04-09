#!/usr/bin/env bash
# =============================================================================
# PlasmidNet EC2 Pipeline: MOBsuite + AMRFinderPlus on 109K NCBI plasmids
#
# Usage:
#   1. Launch EC2 (c6i.4xlarge or larger, 50GB root, Amazon Linux 2023)
#   2. scp this script + data/ncbi_accessions.txt to the instance
#   3. bash ec2_pipeline.sh
#   4. scp results back: results/mob_results.csv, results/amr_results.csv
#   5. Terminate instance
#
# Time: ~2.5h on 16 vCPU, ~45min on 96 vCPU
# =============================================================================
set -euo pipefail

NCPU=$(nproc)
WORKDIR="/home/ec2-user/plasmidnet_pipeline"
ACCESSION_FILE="ncbi_accessions.txt"
BATCH_SIZE=500
FASTA_DIR="$WORKDIR/fasta"
RESULTS_DIR="$WORKDIR/results"
NCBI_API_KEY="${NCBI_API_KEY:-}"  # set for 10 req/sec instead of 3

echo "============================================="
echo "PlasmidNet MOBsuite + AMRFinderPlus Pipeline"
echo "CPUs: $NCPU"
echo "Start: $(date)"
echo "============================================="

mkdir -p "$WORKDIR" "$FASTA_DIR" "$RESULTS_DIR"
cd "$WORKDIR"

# Copy accession file if not already here
if [ ! -f "$ACCESSION_FILE" ]; then
    echo "ERROR: $ACCESSION_FILE not found in $WORKDIR"
    echo "Copy it with: scp data/ncbi_accessions.txt ec2-user@<ip>:$WORKDIR/"
    exit 1
fi

TOTAL_ACC=$(wc -l < "$ACCESSION_FILE")
echo "Accessions to process: $TOTAL_ACC"

# =============================================================================
# Step 1: Install dependencies
# =============================================================================
echo ""
echo ">>> Step 1: Installing dependencies..."

# Install Miniconda if not present
if ! command -v conda &> /dev/null; then
    echo "  Installing Miniconda..."
    curl -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
    conda init bash
    source ~/.bashrc
else
    eval "$(conda shell.bash hook)"
fi

# Create environment with MOBsuite + AMRFinderPlus
if ! conda env list | grep -q plasmidnet; then
    echo "  Creating conda environment..."
    conda create -y -n plasmidnet -c bioconda -c conda-forge \
        mob_suite amrfinderplus parallel curl
fi

conda activate plasmidnet

# Setup AMRFinderPlus database
if [ ! -d "$HOME/amrfinder_db" ]; then
    echo "  Downloading AMRFinderPlus database..."
    amrfinder_update -d "$HOME/amrfinder_db"
fi

echo "  Dependencies ready: $(mob_typer --version 2>&1 | head -1)"
echo "  AMRFinderPlus: $(amrfinder --version 2>&1 | head -1)"

# =============================================================================
# Step 2: Download sequences from NCBI
# =============================================================================
echo ""
echo ">>> Step 2: Downloading sequences from NCBI..."

download_batch() {
    local batch_file="$1"
    local out_file="$2"
    local ids=$(cat "$batch_file" | tr '\n' ',')
    ids=${ids%,}  # remove trailing comma

    local api_param=""
    if [ -n "$NCBI_API_KEY" ]; then
        api_param="&api_key=$NCBI_API_KEY"
    fi

    for attempt in 1 2 3; do
        if curl -sf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ids}&rettype=fasta&retmode=text${api_param}" \
            -o "$out_file" 2>/dev/null; then
            # Verify it's a valid FASTA
            if head -1 "$out_file" | grep -q "^>"; then
                return 0
            fi
        fi
        sleep $((attempt * 2))
    done
    echo "  WARNING: Failed to download batch $(basename $batch_file)" >&2
    return 1
}

export -f download_batch
export NCBI_API_KEY

# Split accessions into batches
split -l $BATCH_SIZE "$ACCESSION_FILE" "$WORKDIR/batch_"
BATCH_COUNT=$(ls "$WORKDIR"/batch_* | wc -l)
echo "  Split into $BATCH_COUNT batches of $BATCH_SIZE"

# Download in parallel (limit to 3 concurrent to respect NCBI rate limits)
CONCURRENT=3
if [ -n "$NCBI_API_KEY" ]; then
    CONCURRENT=8
fi

batch_num=0
for batch_file in "$WORKDIR"/batch_*; do
    batch_num=$((batch_num + 1))
    out_file="$FASTA_DIR/$(basename $batch_file).fasta"
    if [ ! -f "$out_file" ] || [ ! -s "$out_file" ]; then
        echo "$batch_file $out_file"
    fi
done | parallel -j $CONCURRENT --colsep ' ' \
    'download_batch {1} {2} && echo "  Downloaded batch $(basename {1})"' \
    2>&1 | tail -20

# Combine all FASTA files
echo "  Combining FASTA files..."
cat "$FASTA_DIR"/*.fasta > "$WORKDIR/all_plasmids.fasta"
SEQ_COUNT=$(grep -c "^>" "$WORKDIR/all_plasmids.fasta" || echo 0)
FASTA_SIZE=$(du -sh "$WORKDIR/all_plasmids.fasta" | cut -f1)
echo "  Downloaded: $SEQ_COUNT sequences ($FASTA_SIZE)"

# =============================================================================
# Step 3: Run MOBsuite mob_typer
# =============================================================================
echo ""
echo ">>> Step 3: Running MOBsuite mob_typer ($NCPU parallel)..."

# Split combined FASTA into individual files for parallel processing
mkdir -p "$WORKDIR/individual_fasta"

# Use awk to split FASTA by record
awk '/^>/{f="'"$WORKDIR/individual_fasta/"'"++n".fasta"} {print > f}' \
    "$WORKDIR/all_plasmids.fasta"

INDIVIDUAL_COUNT=$(ls "$WORKDIR/individual_fasta/"*.fasta 2>/dev/null | wc -l)
echo "  Split into $INDIVIDUAL_COUNT individual FASTA files"

# Run mob_typer in parallel
run_mob_typer() {
    local fasta="$1"
    local outdir="$2"
    local basename=$(basename "$fasta" .fasta)
    mob_typer -i "$fasta" -o "$outdir/${basename}_mob.txt" 2>/dev/null || true
}
export -f run_mob_typer

mkdir -p "$WORKDIR/mob_output"

ls "$WORKDIR/individual_fasta/"*.fasta | \
    parallel -j "$NCPU" --bar \
    'run_mob_typer {} '"$WORKDIR/mob_output"

# Combine MOBsuite results
echo "  Combining MOBsuite results..."
head -1 "$(ls "$WORKDIR/mob_output/"*_mob.txt 2>/dev/null | head -1)" \
    > "$RESULTS_DIR/mob_results.csv" 2>/dev/null || true

for f in "$WORKDIR/mob_output/"*_mob.txt; do
    tail -n +2 "$f" >> "$RESULTS_DIR/mob_results.csv" 2>/dev/null || true
done

MOB_COUNT=$(wc -l < "$RESULTS_DIR/mob_results.csv")
echo "  MOBsuite results: $((MOB_COUNT - 1)) plasmids typed"

# =============================================================================
# Step 4: Run AMRFinderPlus
# =============================================================================
echo ""
echo ">>> Step 4: Running AMRFinderPlus ($NCPU threads)..."

run_amrfinder() {
    local fasta="$1"
    local outdir="$2"
    local basename=$(basename "$fasta" .fasta)
    amrfinder -n "$fasta" -d "$HOME/amrfinder_db/latest" \
        -o "$outdir/${basename}_amr.tsv" \
        --plus 2>/dev/null || true
}
export -f run_amrfinder

mkdir -p "$WORKDIR/amr_output"

ls "$WORKDIR/individual_fasta/"*.fasta | \
    parallel -j "$NCPU" --bar \
    'run_amrfinder {} '"$WORKDIR/amr_output"

# Combine AMRFinderPlus results
echo "  Combining AMRFinderPlus results..."
head -1 "$(ls "$WORKDIR/amr_output/"*_amr.tsv 2>/dev/null | head -1)" \
    > "$RESULTS_DIR/amr_results.tsv" 2>/dev/null || true

for f in "$WORKDIR/amr_output/"*_amr.tsv; do
    tail -n +2 "$f" >> "$RESULTS_DIR/amr_results.tsv" 2>/dev/null || true
done

AMR_COUNT=$(wc -l < "$RESULTS_DIR/amr_results.tsv")
echo "  AMRFinderPlus results: $((AMR_COUNT - 1)) annotations"

# =============================================================================
# Step 5: Summary
# =============================================================================
echo ""
echo "============================================="
echo "Pipeline Complete!"
echo "End: $(date)"
echo "============================================="
echo ""
echo "Results:"
echo "  MOBsuite:      $RESULTS_DIR/mob_results.csv"
echo "  AMRFinderPlus:  $RESULTS_DIR/amr_results.tsv"
echo ""
echo "Download results:"
echo "  scp ec2-user@\$(curl -s http://169.254.169.254/latest/meta-data/public-ipv4):$RESULTS_DIR/mob_results.csv ."
echo "  scp ec2-user@\$(curl -s http://169.254.169.254/latest/meta-data/public-ipv4):$RESULTS_DIR/amr_results.tsv ."
echo ""
echo "Then on your Mac:"
echo "  python import_mobsuite_results.py"
echo ""
echo "Don't forget to TERMINATE the instance!"
