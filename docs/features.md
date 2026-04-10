# Dashboard Features

PlasmidNet has 10 interactive tabs, each providing a different perspective on the plasmid database.

## Overview

Summary statistics across all 182,362 plasmids:

- **Topology distribution**: 66,511 circular / 6,045 linear
- **Data source**: 59,471 RefSeq / 13,085 INSDC
- **Host kingdom**: 71,949 Bacteria / 607 Archaea
- **Temporal distribution**: Plasmid submissions by year (1982-2024)
- **GC content distribution**: Binned GC% across all plasmids
- **Length distribution**: Plasmid sizes from <5 kb to >500 kb

## Taxonomy

Top 20 host bacterial genera with interactive search.

**Example**: Search for "Klebsiella" to find all 12,213 Klebsiella plasmids in PLSDB.

### Search by genus or species

Enter a genus (e.g., `Escherichia`) or species (e.g., `Escherichia coli`) to get the count from the full database.

## AMR Analysis

Antimicrobial resistance gene frequencies across the database.

### Drug Class Distribution
- Top 20 AMR drug classes (excluding heavy metals, shown separately)
- Includes: Beta-lactam, Aminoglycoside, Sulfonamide, Tetracycline, Quinolone, Trimethoprim, Phenicol, QAC

### Individual Gene Frequency
- Top 30 AMR genes coloured by their real drug class from the database
- **Key genes**: sul1 (5,901), qacEdelta1 (5,446), blaTEM-1 (5,250), sul2 (4,379)

## Plasmid Viewer

Interactive circular plasmid map for any NCBI accession.

### How to use

1. Enter an NCBI accession (e.g., `MH595533`) or `NZ_CP031107.1`
2. Click "Visualize"
3. Scroll to zoom, hover for gene details

### Ring layout (outside-in)

| Ring | Content | Color |
|------|---------|-------|
| Outer | Forward-strand CDS | Blue (AMR in red) |
| Middle | Reverse-strand CDS | Green (AMR in red) |
| Inner features | AMR / MOB / BGC from PLSDB | Red / Amber / Purple |
| Center | GC content deviation | Teal (above mean) / Red (below) |

### Gene classification

Genes are automatically classified:
- **AMR**: bla*, aph*, aac*, tet*, sul*, dfr*, qnr*, mcr*, fos*
- **MOB**: mob*, tra*, trb*, virB*, relaxase
- **TA**: relE, higA, ccdA/B, vapB, pemK
- **VIR**: vir*, iro*, iuc*, hly*, tsh*
- **QAC**: qacE, qacA, smr, emrE
- **Metal**: mer*, ars*, pco*, sil*, ter*

### Example: MH595533 (KPC plasmid)

```
Accession: MH595533
Length: 14,873 bp | circular | GC 54.6%
Organism: Klebsiella pneumoniae

Detected:
  22 CDS genes
  2 AMR: aphA(3')-Vla, blaKPC-2
  2 MOB: mobC, mobA
  3 TA: relE, higA, higB
```

## Inc/MOB/IS

### Incompatibility Groups
Top 20 replicon types. IncFIB and IncFII dominate, followed by IncFIA, IncI1, ColRNAI.

### Predicted Mobility
- **Conjugative**: 22,588 (31%) — can self-transfer
- **Mobilizable**: 21,152 (29%) — needs helper plasmid
- **Non-mobilizable**: 28,816 (40%)

### Relaxase and MPF Types
Relaxase classification (MOBF, MOBP, MOBQ, MOBH, MOBC, MOBV) and Mating Pair Formation types (MPF_F, MPF_T, MPF_I).

### IS Elements & Transposon Families
Distribution of 438,616 IS/transposon annotations:
- **IS26**: 48,228 annotations (13,298 plasmids) — dominant in AMR gene cassettes
- **IS3**: 21,125 | **IS5**: 14,932 | **IS110**: 10,905 | **IS6**: 9,960
- **Tn3**: 9,361 | **IS1**: 8,136 | **IS30**: 6,961

## Correlations

Cross-feature co-occurrence analysis using the **entire 72,556-plasmid PLSDB database**.

### 1. AMR Drug Class vs Inc Groups & Mobility
Heatmap showing which Inc groups carry which resistance classes. IncFIB dominates beta-lactam and aminoglycoside resistance.

### 2. Heavy Metal Resistance
Mercury, arsenic, copper, silver, and tellurium resistance genes. IncHI2 is enriched for tellurium resistance.

### 3. Phage & Mobile Elements
Distribution of phage and transposase elements per plasmid, correlated with mobility type.

### 4. pMLST Distribution
Plasmid MLST scheme frequency and allele distribution.

### 5. Virulence Factors (6,145 plasmids)
Top virulence genes: clpK, hsp20, iroD, iroN, iucC. 58% on conjugative plasmids.

### 6. Toxin-Antitoxin Systems (40,759 plasmids)
Top TA genes from PGAP: vapB (6,569), ccdA (4,077), ccdB (4,036), pemI (3,554).

### 7. QAC Resistance (7,023 plasmids)
Quaternary ammonium compound resistance. 63% on conjugative plasmids.

## Geography

57,751 geolocated plasmids across 100+ countries.

### 1. Global Map
Scatter map coloured by predicted mobility. Hover for accession, country, and Inc group.

### 2. Country Comparison
- **Mobility by country**: China leads (9,412 plasmids), followed by USA (8,821)
- **Inc groups by country**: Heatmap showing which Inc groups dominate in each country

### 3. Temporal Spread Animations
Animated timeline (2000-2024) showing cumulative plasmid reports per country. Press Play to watch the global spread. Separate animations for specific Inc groups and mobility types.

### 4. Plasmid Host & Source
Distribution by host category:
- **Human**: 11,464 | **Animal**: 5,925 | **Soil**: 1,576
- **Water**: 1,332 | **Food**: 1,059 | **Environment**: 959

Heatmaps for Inc groups and AMR drug classes by host source.

## Analytics

Statistical analyses for confound decomposition.

### 1. Matched Comparison
Conjugative % by species x country. Controls for species composition bias.

### 2. Rarefaction Analysis
Bootstrap subsampling to test whether mobility differences between countries are real or artifacts of sample size.

### 3. XGBoost + SHAP
XGBoost classifier (82.6% accuracy) predicting mobility from country, genus, host source, year, length, GC%, and Inc groups. SHAP values show feature importance.

### 4. Temporal Trends
AMR drug class and Inc group prevalence over time (2010-2024).

### 5. Simpson's Paradox Detector
Flags Inc groups where the overall mobility trend reverses when stratified by species.

### 6. AMR Co-occurrence Network
Heatmap of AMR genes co-occurring on the same plasmid. Reveals multi-drug resistance cassettes.

### 7. Integron & Gene Cassette Analysis
Class 1 integrons (5,371 plasmids with qacEdelta1 + sul1). Gene cassettes within 5 kb of qacEdelta1: sul1, aadA2, dfrA12, arr-3, catB3, blaOXA-1.

### 8. Co-mobilization Analysis
45% of mobilizable plasmids share a host with a conjugative plasmid. Relaxase x T4SS compatibility heatmap. ML predictor (77% accuracy) for co-mobilization.

## Compare

Compare up to 10 plasmids with BLASTn alignment and pyCirclize visualization.

### How to use

1. Enter NCBI accessions (one per line) or upload FASTA files
2. Click "Compare"
3. View circular comparison plot, alignment table, and CDS annotations

### Example

```
MH595534
MH595533
NC_011102.1
```

### Output
- **Circular plot**: One ring per plasmid with CDS genes and GC content
- **Alignment links**: Colour-coded by identity (green >99%, blue 95-99%)
- **CDS tables**: Gene annotations for each plasmid
- **Downloads**: FASTA and GenBank format for each plasmid

## Seq Analysis

Scan any DNA sequence for restriction sites, promoters, RBS, codon bias, IS elements, and engineering indicators.

### How to use

1. Enter an NCBI accession, upload a FASTA file, or paste raw DNA
2. Click "Fetch & Analyze" or "Analyze Sequence"
3. Review the engineering score and detailed findings

### Features detected

| Feature | Method |
|---------|--------|
| Restriction sites | 25 common enzymes + reverse complements |
| RE hotspots | 500 bp windows with 3+ sites |
| Cryptic promoters | -35 + spacer + -10 (sigma-70) |
| RBS (Shine-Dalgarno) | AAGGAG, AGGAGG + downstream ATG |
| Vector signatures | pUC, pBR322, ColE1, T7/T3/SP6, CMV, lacZ, f1 |
| IS elements | IS1, IS26, IS903, IS10, IS3, IS5, ISEcp1, Tn3, Tn21 |
| Direct repeats | 4-12 bp TSDs flanking insertions |
| NTEKPC detection | Tn4401 vs NTEKPC variants around blaKPC |
| Codon usage bias | Per-window CAI vs E. coli K-12 |
| K-mer naturalness | 4-mer entropy + palindrome fraction |

### Engineering Score (0-100)

| Score | Classification |
|-------|---------------|
| 0-19 | Natural |
| 20-39 | Natural (resistance platform) |
| 40-59 | Ambiguous |
| 60-79 | Likely engineered |
| 80-100 | Engineered |

### Mobilization Assessment

Classifies plasmids by transfer capability:
- **Self-transmissible**: Has relaxase + T4SS
- **Mobilizable**: Has relaxase, no T4SS
- **Retro-mobilizable**: Only oriT, needs all machinery in trans
- **Non-mobilizable**: No transfer genes

### Example: MH595533

```
Score: 45/100 (Ambiguous)
Classification: Natural resistance platform

NTEKPC: blaKPC detected at position 12,829
Direct repeats: 7 TSDs found
Mobilization: Mobilizable (MOB genes) — has mobA + mobC, no T4SS
Vector signatures: pBR322_ori, ColE1_ori (shared natural/synthetic)
IS elements: 0 (ISEcp1-like element not in current signature set)
```
