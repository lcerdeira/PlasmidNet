# Dashboard Features

PlasmidNet has 11 interactive tabs providing comprehensive plasmid analytics across **182,362 plasmids** (72,556 PLSDB + 109,806 NCBI) with **600,896 AMR annotations**.

## Overview

Hero banner with total plasmid count, source breakdown (PLSDB + NCBI), and key statistics. Collapsible "About PlasmidNet" section with project motivation.

Charts: topology, source, kingdom distributions, temporal trends (1982-2024), GC content, and plasmid length distributions.

## Taxonomy

Top 20 host bacterial genera with interactive search by genus or species. Top genera: Escherichia (14,804), Klebsiella (12,213), Salmonella (3,524), Enterococcus (3,231), Staphylococcus (3,180).

## AMR

- **Drug Class Distribution**: Top 20 AMR drug classes from the full database
- **Individual Gene Frequency**: Top 30 AMR genes coloured by real drug class (sul1, qacEdelta1, blaTEM-1, sul2, tet(A), aph(6)-Id)

## Viewer

Interactive circular plasmid map (Plotly Barpolar) for any NCBI accession with zoom and hover tooltips.

Genes auto-classified: AMR (red), MOB (amber), TA (orange), VIR (dark red), QAC (teal), Metal (purple). GC content deviation ring. Works for PLSDB and any NCBI GenBank accession (automatic fallback).

## Inc/MOB/IS

### Incompatibility Groups & Mobility
- Top 20 replicon types (IncFIB, IncFII dominant)
- Predicted mobility: conjugative 31%, mobilizable 29%, non-mobilizable 40%
- Relaxase types (MOBF, MOBP, MOBQ, MOBH, MOBC, MOBV)
- MPF types (MPF_F, MPF_T, MPF_I)

### IS Elements & Transposon Families (438,616 annotations)
- IS26 (48,228), IS3 (21,125), IS5 (14,932), IS110 (10,905), Tn3 (9,361)
- IS families by mobility type
- IS families by Inc group (heatmap)

## Correlations (n = 72,556)

1. **AMR Drug Class vs Inc Groups & Mobility**
2. **Heavy Metal Resistance** (mercury, arsenic, copper, silver, tellurium)
3. **Phage & Mobile Elements**
4. **pMLST Distribution**
5. **Virulence Factors** (6,145 plasmids)
6. **Toxin-Antitoxin Systems** (40,759 plasmids from PGAP)
7. **QAC Resistance** (16,948 plasmids)

## Geography (57,751 geolocated plasmids)

1. **Global Plasmid Map** by mobility (scatter map with hover)
2. **Country Comparison** — mobility and Inc distribution by top 20 countries
3. **Temporal Spread Animations** — animated timeline 2000-2024 for all plasmids, specific Inc groups, and mobility types (dropdown-driven)
4. **Plasmid Host & Source** — Human (11,464), Animal (5,925), Soil (1,576), Water (1,332), Food (1,059), Environment (959). Heatmaps for Inc groups and AMR by host.

## Analytics

### UMAP Clustering
15,000 plasmids in 2D feature space coloured by mobility and genus. Reveals natural plasmid ecotypes.

### Matched Comparison
Conjugative % by species x country. Controls for species composition bias. Key finding: Escherichia in China 54% conjugative vs 35% in Germany.

### Rarefaction Analysis
Bootstrap subsampling (50x) to test sample size effects on mobility ratios.

### XGBoost + SHAP (82.6% accuracy)
Pre-computed at startup. Feature importance for predicting mobility.

### Temporal Trends (2010-2024)
AMR drug class and Inc group prevalence over time.

### Simpson's Paradox Detector
Flags Inc groups with species-level trend reversals (>20pp spread).

### AMR Co-occurrence Network
Heatmap of resistance genes co-occurring on the same plasmid.

### Integron & Gene Cassette Analysis (5,371 plasmids)
Class 1 integrons (qacEdelta1 + sul1). Gene cassettes within 5kb: sul1, aadA2, dfrA12, arr-3, catB3, blaOXA-1. 74% on conjugative plasmids.

### Co-mobilization (45% co-location rate)
Conjugative x mobilizable Inc group pairs. Relaxase x T4SS compatibility heatmap. ML predictor (77% accuracy) identifying relaxase type as strongest predictor.

### Retro-mobilization & HGT Routes
Transfer mechanisms for 28,816 non-mobilizable plasmids:
- **Retro-mobilization**: 9,555 (33%) share host with conjugative plasmid
- **Mobilizable relay**: 10,736 (37%)
- **Transduction**: 8,585 (30%) small enough for phage capsids
- **2,226 non-mob with AMR + conjugative partner** = retro-mobilization risk

### blaKPC Transposon Context (4,800+ plasmids)
Genes within 5kb of blaKPC (blaTEM, blaCTX-M-65, sul1). NTEKPC elements characterised as true transposons. blaKPC by Inc group (IncFIB/FII, IncN) and mobility (59% conjugative).

### Integron ML & Transposon-AMR Correlations
ML predicting integron carriage. IS family x AMR drug class co-occurrence heatmap.

## Compare

Compare up to 10 plasmids with BLASTn + pyCirclize:
- Enter NCBI accessions or upload FASTA files
- CDS gene annotations on circular plot (AMR red, MOB amber)
- Alignment links colour-coded by identity
- CDS annotation tables per plasmid
- Download as FASTA or GenBank format

## Seq Analysis

Scan any DNA sequence (NCBI accession, file upload, or paste/edit):

| Feature | Method |
|---------|--------|
| Restriction sites | 25 common enzymes + reverse complements |
| RE hotspots | 500bp windows with 3+ sites |
| Cryptic promoters | -35 + spacer + -10 (sigma-70) |
| RBS (Shine-Dalgarno) | AAGGAG, AGGAGG + downstream ATG |
| Vector signatures | pUC, pBR322, ColE1, T7/T3/SP6, CMV, lacZ, f1 |
| IS elements | IS1, IS26, IS903, IS10, IS3, IS5, ISEcp1, Tn3, Tn21 |
| Direct repeats / TSDs | 4-12bp flanking insertions |
| NTEKPC detection | Tn4401 vs NTEKPC variants around blaKPC |
| Codon usage bias | Per-window CAI vs E. coli K-12 |
| K-mer naturalness | 4-mer entropy + palindrome fraction |
| Engineering score | 0-100 with 5-tier classification |
| Mobilization assessment | Retro-mobilizable / mobilizable / conjugative |

## API

Interactive API Explorer tab with dropdown + query input. Also accessible programmatically:

| Endpoint | Description |
|----------|-------------|
| `GET /api/` | Documentation |
| `GET /api/plasmid/<acc>` | Full plasmid details |
| `GET /api/search?q=<query>` | Search plasmids |
| `GET /api/amr?gene=<gene>` | AMR gene search |
| `GET /api/amr/classes` | Drug class counts |
| `GET /api/typing?inc=<inc>` | Inc group search |
| `GET /api/typing/mobility` | Mobility distribution |
| `GET /api/taxonomy?genus=<genus>` | Taxonomy search |
| `GET /api/stats` | Database summary |
