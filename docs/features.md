# Dashboard Features

PlasmidNet has 10 interactive tabs providing comprehensive plasmid analytics.

## Overview

Summary statistics across **182,362 plasmids** (72,556 PLSDB + 109,806 NCBI):

- **Topology**: 66,511 circular / 6,045 linear
- **Source**: 59,471 RefSeq / 13,085 INSDC
- **Host kingdom**: 71,949 Bacteria / 607 Archaea
- **AMR annotations**: 600,896 total (251K PLSDB + 350K NCBI)
- **Temporal distribution**: 1982-2024
- **GC content** and **length distributions**

## Taxonomy

Top 20 host bacterial genera with interactive search by genus or species.

**Top genera**: Escherichia (14,804), Klebsiella (12,213), Salmonella (3,524), Enterococcus (3,231), Staphylococcus (3,180)

## AMR Analysis

- **Drug Class Distribution**: Top 20 AMR drug classes (Beta-lactam, Aminoglycoside, Sulfonamide, etc.)
- **Individual Gene Frequency**: Top 30 AMR genes coloured by real drug class from database
- **Key genes**: sul1 (5,901), qacEdelta1 (5,446), blaTEM-1 (5,250)

## Plasmid Viewer

Interactive circular plasmid map for any NCBI accession with zoom and hover tooltips.

### Gene classification (automatic)
- **AMR**: bla*, aph*, aac*, tet*, sul*, dfr*, qnr*, mcr*, fos*
- **MOB**: mob*, tra*, trb*, virB*, relaxase
- **TA**: relE, higA, ccdA/B, vapB, pemK
- **VIR**: vir*, iro*, iuc*, hly*, tsh*
- **QAC**: qacE, qacA, smr, emrE
- **Metal**: mer*, ars*, pco*, sil*, ter*

## Inc/MOB/IS

### Incompatibility Groups
Top 20 replicon types. IncFIB and IncFII dominate.

### Predicted Mobility
- Conjugative: 22,588 (31%)
- Mobilizable: 21,152 (29%)
- Non-mobilizable: 28,816 (40%)

### Relaxase and MPF Types
MOBF, MOBP, MOBQ, MOBH, MOBC, MOBV and MPF_F, MPF_T, MPF_I.

### IS Elements & Transposon Families (438,616 annotations)
- IS26 (48,228), IS3 (21,125), IS5 (14,932), IS110 (10,905), IS6 (9,960)
- Tn3 (9,361), IS1 (8,136), IS30 (6,961)
- Heatmap: IS families by Inc group

## Correlations (n = 72,556)

1. **AMR Drug Class vs Inc Groups & Mobility** — heatmap and grouped bar
2. **Heavy Metal Resistance** — mercury, arsenic, copper, silver, tellurium
3. **Phage & Mobile Elements** — distribution and mobility correlation
4. **pMLST Distribution** — scheme frequency, alleles, mobility
5. **Virulence Factors** (6,145 plasmids) — clpK, hsp20, iroD, iucC
6. **Toxin-Antitoxin** (40,759 plasmids) — vapB, ccdA/B, pemI/K from PGAP
7. **QAC Resistance** (16,948 plasmids) — quaternary ammonium compounds

## Geography (57,751 geolocated plasmids)

1. **Global Plasmid Map** — scatter map by mobility
2. **Country Comparison** — mobility and Inc distribution by top 20 countries
3. **Temporal Spread Animations** — animated timeline 2000-2024 for all plasmids, specific Inc groups, and mobility types
4. **Plasmid Host & Source** — Human (11,464), Animal (5,925), Soil (1,576), Water (1,332), Food (1,059), Environment (959)

## Analytics

### UMAP Clustering
15,000 plasmids embedded in 2D feature space by Inc group, genus, host, length, GC. Reveals natural plasmid ecotypes.

### Matched Comparison
Conjugative % by species x country (controls for species composition).

### Rarefaction Analysis
Bootstrap subsampling to test sample size effects.

### XGBoost + SHAP (82.6% accuracy)
Feature importance for predicting mobility. SHAP values show each feature's contribution.

### Temporal Trends (2010-2024)
AMR drug class and Inc group prevalence over time.

### Simpson's Paradox Detector
Flags Inc groups with species-level trend reversals.

### AMR Co-occurrence Network
Which resistance genes travel together on the same plasmid.

### Integron & Gene Cassette Analysis (5,371 plasmids)
Class 1 integrons: qacEdelta1 + sul1 + aadA gene cassettes. 74% on conjugative plasmids.

### Co-mobilization (45% co-location rate)
Relaxase x T4SS compatibility heatmap. ML predictor (77% accuracy).

## Compare

Compare up to 10 plasmids with BLASTn alignment and pyCirclize visualization.
- Upload FASTA files or enter NCBI accessions
- CDS gene annotations on circular plot
- Alignment links colour-coded by identity
- Download as FASTA or GenBank

## Seq Analysis

Scan any DNA sequence for:
- Restriction sites (25 enzymes) and RE hotspots
- Cryptic promoters (-35 + -10 sigma-70)
- Ribosome Binding Sites (Shine-Dalgarno)
- Vector backbone signatures (pUC, pBR322, ColE1, T7, CMV)
- IS elements and transposons
- Direct repeats / target site duplications
- NTEKPC transposon detection (Tn4401 vs NTEKPC variants)
- Codon usage bias and k-mer naturalness
- Engineering score (0-100) with 5-tier classification
- Mobilization assessment (retro-mobilization detection)

## REST API

Programmatic access at `/api/`:
- `/api/plasmid/<acc>` — full plasmid details
- `/api/search?q=<query>` — search by description/organism
- `/api/amr?gene=<gene>` — AMR gene search
- `/api/typing?inc=<inc>` — Inc group search
- `/api/stats` — database summary
