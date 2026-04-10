# Plasmid Comparison

Compare up to 10 plasmids with BLASTn alignment and pyCirclize circular visualization.

## How to Use

1. **Enter accessions** (one per line) in the text area:
   ```
   MH595534
   MH595533
   NC_011102.1
   ```
2. **Upload FASTA files** by dragging & dropping or clicking the upload area
3. **Mix both**: Combine accessions with uploaded files (up to 10 total)
4. Click **Compare**

## Output

### Circular Comparison Plot

Generated with [pyCirclize](https://github.com/moshi4/pyCirclize):
- One sector per plasmid with distinct colours
- **Outer ring**: CDS genes (blue = CDS forward, green = reverse, red = AMR, amber = MOB)
- **Inner ring**: GC content deviation from mean
- **Links between sectors**: BLASTn alignment regions, colour-coded by identity

### Alignment Table

BLASTn HSPs (High-Scoring Segment Pairs):
- Query and subject accessions
- Percent identity
- Alignment length
- Query and subject coordinates

### CDS Annotation Tables

For each plasmid with GenBank annotations:
- Gene name, product description
- Start/end coordinates, strand

### Downloads

- **FASTA**: Raw sequence for each plasmid
- **GenBank**: Annotated GenBank format with CDS features

## Requirements

- **BLASTn** must be installed for alignment (`ncbi-blast+` package)
- On Heroku, installed via the `Aptfile` buildpack
- Locally: `brew install blast` (macOS) or `apt install ncbi-blast+` (Ubuntu)

## Example Comparison

Comparing four plasmids from the same study:

```
MH595534  (9,548 bp, 11 CDS)
MH595533  (14,873 bp, 22 CDS — carries blaKPC-2)
MG786907  (3,812 bp, no annotations)
MH000708  (6,367 bp, no annotations)
```

The comparison reveals:
- MH595533 and MH595534 share regions of high identity (>95%)
- MG786907 and MH000708 lack CDS annotations but their sequences are included in BLASTn
- Gene names (blaKPC-2, mobA, repA, etc.) are labelled on the outer ring
