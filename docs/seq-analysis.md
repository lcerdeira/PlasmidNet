# Sequence Analysis

The Seq Analysis tab scans any DNA sequence for functional elements, cloning artifacts, and engineering indicators.

## Input Methods

1. **NCBI accession**: Enter an accession (e.g., `MH595533`) and click "Fetch & Analyze"
2. **File upload**: Drag & drop a FASTA or GenBank file
3. **Paste/Edit**: Type or paste directly into the textarea, modify, then click "Analyze Sequence"

## Analysis Pipeline

### 1. Restriction Site Mapping

Scans for 25 common Type II restriction enzymes including:
- **6-cutters**: EcoRI, BamHI, HindIII, SalI, XhoI, NdeI, NcoI, XbaI, PstI, SacI, KpnI, NotI
- **Type IIS (Golden Gate)**: BsaI, BbsI, BsmBI, SapI
- Both forward and reverse complement orientations

**Hotspot detection**: 500 bp sliding windows with 3+ restriction sites indicate potential cloning junctions.

### 2. Promoter Prediction

Sigma-70 promoter consensus search:
- **-35 box**: TTGAC[AT]
- **Spacer**: 14-20 bp (optimal 17 bp)
- **-10 Pribnow box**: TA[AT]AAT

Scored 0-100 based on consensus match and spacer length.

### 3. Ribosome Binding Sites

Shine-Dalgarno sequence detection:
- **Strong**: AAGGAG, AGGAGG
- **Moderate**: AGGAG
- **Weak**: GAGG

Checks for downstream ATG start codon within 5-13 bp.

### 4. Vector Backbone Signatures

Detects sequences from common cloning vectors:

| Signature | Source | Score impact |
|-----------|--------|-------------|
| T7/T3/SP6 promoter | Synthetic only | +15 (strong) |
| lacZ alpha | Synthetic only | +15 (strong) |
| f1 origin | Synthetic only | +15 (strong) |
| CMV promoter | Synthetic only | +15 (strong) |
| ColE1/pBR322/pUC origin | Shared natural/synthetic | +3 (weak) |

### 5. IS Element & Transposon Detection

Scans for terminal inverted repeats of common IS families:
- IS1, IS26, IS903, IS10, IS3, IS5, ISEcp1
- Tn3, Tn21, Tn1721

IS elements near RE hotspots **reduce** the engineering score (natural transposon boundaries, not cloning junctions).

### 6. Direct Repeats / Target Site Duplications

Finds 4-12 bp direct repeats separated by 500-5000 bp. These are created when a transposon inserts into DNA and are hallmarks of natural mobile element activity.

### 7. NTEKPC Transposon Detection

For plasmids carrying blaKPC:
- **Tn4401**: Classic ISKpn7-blaKPC-ISKpn6 structure
- **NTEKPC-IId**: ISKpn27-blaKPC-tnpA structure
- **NTEKPC variants**: Other non-Tn4401 contexts

### 8. Codon Usage Analysis

Per-window Codon Adaptation Index (CAI) relative to E. coli K-12. High CAI regions (>0.7) may indicate codon optimization. Low CAI regions (<0.2) indicate unusual codon usage.

### 9. K-mer Naturalness

4-mer frequency analysis:
- **Entropy ratio**: Natural DNA has higher entropy than engineered
- **Palindrome fraction**: High palindrome content correlates with restriction site density
- **Perfect repeats**: Long identical repeats (>20 bp) are common in synthetic constructs

### 10. Mobilization Assessment

Classifies plasmid transfer capability from CDS annotations:

| Category | Criteria | Example |
|----------|----------|---------|
| Self-transmissible | Relaxase + T4SS | Large IncF conjugative plasmids |
| Mobilizable | Relaxase, no T4SS | Many small ColE1 plasmids |
| Mobilizable (MOB) | MOB genes, no T4SS | MH595533 (mobA + mobC) |
| Retro-mobilizable | Only oriT | Can only transfer with helper |
| Non-mobilizable | No transfer genes | Cryptic plasmids |

## Engineering Score

Composite score (0-100) calibrated against the PLSDB natural plasmid baseline:

- **Engineering scars** (paired BsaI sites, MCS): +20-40 per scar
- **Synthetic-specific vector signatures**: +15 each
- **Shared natural/synthetic signatures**: +3 each
- **RE density above PLSDB 95th percentile**: +5-15
- **Unexplained RE hotspots** (not near IS elements): +5-15
- **IS elements present**: -5 to -15 (natural indicator)
- **K-mer deviation**: +0-20

### Interpretation

| Range | Classification | Typical examples |
|-------|---------------|-----------------|
| 0-19 | Natural | Wild-type resistance plasmids |
| 20-39 | Natural (resistance platform) | Multi-drug resistance plasmids with natural origins |
| 40-59 | Ambiguous | Plasmids with ColE1-type origins (shared with lab vectors) |
| 60-79 | Likely engineered | Plasmids with synthetic promoters or MCS |
| 80-100 | Engineered | Lab vectors (pUC19, pET, etc.) |

## Case Study: MH595533 (pKPN535a)

This 14,873 bp IncQ1 plasmid carries blaKPC-2 and was isolated from Klebsiella pneumoniae.

**Score: 45/100 (Ambiguous)**

Why ambiguous, not engineered?
- ColE1/pBR322 origins are **shared** between natural IncQ plasmids and lab vectors (+3 each, not +15)
- No synthetic-specific signatures (no T7 promoter, no lacZ, no f1 ori)
- No IS elements detected (the ISKpn27 flanking blaKPC uses a different signature)
- Direct repeats found (7 TSDs) — evidence of natural transposon activity
- Mobilization: mobA + mobC present, classified as "Mobilizable (MOB genes)"
