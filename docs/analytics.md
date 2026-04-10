# Statistical Analytics

The Analytics tab provides confound decomposition to answer: **are regional differences in plasmid mobility biological or sampling bias?**

## 1. Matched Comparison

**Question**: Do mobility differences persist when comparing the same species across countries?

**Method**: Heatmap of conjugative % for each (species, country) pair. Green = more conjugative, red = less.

**Key finding**: *Escherichia* in China is 54% conjugative vs 35% in Germany — the effect persists within the same species, suggesting a real biological or sampling-context difference.

## 2. Rarefaction Analysis

**Question**: Are mobility differences artifacts of sample size?

**Method**: Each country subsampled to equal n (50 to 2000), bootstrapped 50 times. Error bars show standard deviation.

**Interpretation**: If curves for different countries converge at the same value → differences are noise. If they stabilize at different values → differences are real.

## 3. XGBoost + SHAP

**Question**: What features best predict mobility?

**Method**: XGBoost classifier trained on 57,580 plasmids. Features: country, genus, host source, year, length, GC%, Inc groups.

**Results**: 82.6% accuracy. SHAP values show feature importance.

If "Country" ranks high after controlling for species, host source, and year → regional differences are independent.

## 4. Temporal Trends

**Question**: Are AMR genes and Inc groups increasing or declining over time?

**Method**: Line charts showing % prevalence per year (2010-2024), normalized by total plasmids submitted that year.

## 5. Simpson's Paradox Detector

**Question**: Does the overall mobility trend for an Inc group reverse when stratified by species?

**Method**: For each Inc group with >50 plasmids, compare overall conjugative % to per-species conjugative %. Flag cases with >20 percentage point divergence.

**Finding**: No true paradoxes detected — regional mobility differences are consistent across species.

## 6. AMR Co-occurrence Network

**Question**: Which resistance genes travel together?

**Method**: Heatmap of pairwise gene co-occurrence on the same plasmid (minimum 100 co-occurrences). Reveals multi-drug resistance cassettes.

## 7. Integron & Gene Cassette Analysis

**Question**: What gene cassettes are carried by class 1 integrons?

**Method**: Identify plasmids with qacEdelta1 + sul1 (class 1 integron 3' conserved segment). Map all genes within 5 kb of qacEdelta1.

**Finding**: 5,371 plasmids carry class 1 integrons. Top cassettes: sul1, aadA2, dfrA12, arr-3, catB3, blaOXA-1. 74% on conjugative plasmids.

## 8. Co-mobilization Analysis

**Question**: Which mobilizable plasmids are likely co-transferred with conjugative plasmids?

**Findings**:
- **45%** of mobilizable plasmids share a host with a conjugative plasmid
- IncFII conjugative plasmids most frequently carry ColRNAI mobilizable plasmids
- ML predictor (77% accuracy) identifies relaxase type as the strongest predictor of co-mobilization

### Relaxase x T4SS Compatibility

Heatmap of observed co-location vs literature-known compatibility:
- MOBF → MPF_F (known)
- MOBP → MPF_T, MPF_F (known)
- MOBC → MPF_F (885 co-occurrences, **not in standard rules** — potential novel pathway)
