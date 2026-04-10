# PlasmidNet: an interactive platform for comprehensive plasmid analytics integrating 182,362 bacterial plasmids with machine learning-driven insights

Louise Cerdeira^1*^

^1^ London School of Hygiene & Tropical Medicine, London, United Kingdom

*Corresponding author: louise.cerdeira@lshtm.ac.uk

## Abstract

Plasmids are central drivers of antimicrobial resistance (AMR) dissemination, yet systematic tools for exploring their diversity, mobility, and geographic spread remain limited. Here we present PlasmidNet, an open-source interactive dashboard that integrates 182,362 complete bacterial plasmids from PLSDB and NCBI Nucleotide into a unified analytical platform. PlasmidNet provides real-time querying of 600,896 AMR annotations, 613,902 PGAP feature annotations, and geographic metadata for 57,751 geolocated plasmids across 100+ countries. The platform incorporates machine learning models for mobility prediction (82.6% accuracy), co-mobilization assessment (77% accuracy), and integron carriage prediction, alongside statistical analyses for confound decomposition including rarefaction, matched comparisons, and Simpson's paradox detection. Novel features include retro-mobilization risk assessment for non-mobilizable plasmids, NTEKPC transposon context analysis, sequence-level engineering detection with IS element-aware scoring, and a REST API for programmatic access. PlasmidNet is freely available at https://plasmidnet-7302570f4818.herokuapp.com/ with source code at https://github.com/lcerdeira/plasmidnet.

**Keywords**: plasmid, antimicrobial resistance, horizontal gene transfer, conjugation, mobilization, database, dashboard, machine learning

## Introduction

Antimicrobial resistance (AMR) represents one of the most pressing global health challenges, with plasmid-mediated horizontal gene transfer (HGT) serving as the predominant mechanism for resistance gene dissemination among bacterial populations (1, 2). The World Health Organization has identified carbapenem-resistant Enterobacteriaceae as critical priority pathogens, with plasmid-borne carbapenemases such as blaKPC, blaNDM, and blaOXA-48 driving treatment failures worldwide (3). Understanding the epidemiology and dynamics of resistance plasmids requires comprehensive analytical tools that can integrate diverse data types across genomic, geographic, and temporal dimensions.

Several databases have been developed to catalogue bacterial plasmids. PLSDB, the most comprehensive curated plasmid database, contains 72,556 entries with MOB typing, AMR annotations, and taxonomic classification (4, 5). PlasmidFinder enables replicon typing through sequence-based identification of incompatibility groups (6), while MOBsuite provides computational classification of plasmid mobility based on relaxase, mating pair formation (MPF), and origin of transfer (oriT) detection (7). The Comprehensive Antibiotic Resistance Database (CARD) and AMRFinderPlus provide standardised resistance gene annotation (8, 9). Despite the availability of these resources, researchers face significant challenges in synthesising information across platforms, particularly when investigating questions that span AMR prevalence, geographic distribution, plasmid mobility, and host specificity simultaneously.

The role of mobile genetic elements — insertion sequences (IS), transposons, and integrons — in shaping plasmid architecture and facilitating AMR gene mobilisation has received increasing attention (10, 11). Class 1 integrons, characterised by the intI1 integrase and the 3' conserved segment containing qacEdelta1 and sul1, serve as platforms for gene cassette acquisition and are strongly associated with multi-drug resistance phenotypes (12). The interplay between IS elements, transposons such as Tn4401 and the non-Tn4401 elements (NTEKPC) carrying blaKPC, and their host plasmid backbones determines the dissemination dynamics of critical resistance genes (13, 14).

Recent advances in machine learning have demonstrated the potential for predictive modelling in plasmid biology. Gradient boosting methods and random forests have been applied to predict plasmid host range, mobility classification, and AMR gene carriage (15, 16). SHAP (SHapley Additive exPlanations) values provide interpretable feature importance estimates that can distinguish biological signals from sampling artifacts (17). Dimensionality reduction techniques such as UMAP (Uniform Manifold Approximation and Projection) enable visualisation of high-dimensional plasmid feature spaces, revealing natural clustering patterns that correspond to functional plasmid ecotypes (18).

A critical but underexplored aspect of plasmid epidemiology is the transfer potential of non-mobilizable plasmids. While approximately 40% of sequenced plasmids lack autonomous transfer machinery, multiple HGT mechanisms — including retro-mobilization (reverse transfer via conjugative T4SS), conduction (co-integration via shared IS elements), and phage-mediated transduction — can facilitate their dissemination (19, 20). Understanding which non-mobilizable plasmids are at risk of mobilisation, and by which mechanisms, has direct implications for predicting AMR spread in clinical and environmental settings.

Here we present PlasmidNet, an integrated analytical platform that addresses these challenges by combining the largest unified plasmid dataset (182,362 entries) with interactive visualisation, machine learning-driven analytics, and a programmatic REST API. PlasmidNet enables researchers to explore plasmid diversity, AMR patterns, geographic distribution, and transfer mechanisms through a single, freely accessible web interface.

## Materials and Methods

### Data integration and database construction

PlasmidNet integrates data from two primary sources. The PLSDB 2025 database (version 2024_05_31_v2) was obtained from Figshare (4), comprising 72,556 plasmid entries with associated metadata including nucleotide accession, sequence length, GC content, topology, MOB typing results (replicon type, relaxase type, MPF type, oriT type, predicted mobility), pMLST classification, taxonomy, biosample location, and associated PubMed identifiers. AMR annotations (251,138 entries) were derived from AMRFinderPlus analysis provided by PLSDB, covering gene symbols, drug classes, antimicrobial agents, and genomic coordinates. PGAP protein annotations were filtered to extract 613,902 entries corresponding to toxin-antitoxin (TA) systems, phage-related elements, and transposase/IS element families. Geographic metadata was obtained from the biosample table, with country names normalised using latitude/longitude-based reverse geocoding for entries with ambiguous location strings.

An additional 109,806 complete plasmid sequences not present in PLSDB were harvested from NCBI Nucleotide using the Entrez E-utilities API. Duplicate detection was performed by matching both full accession versions and base accessions (without version suffix) against the existing PLSDB entries. AMR annotations for these additional plasmids were generated using AMRFinderPlus v3.12.8 (9) with the 2024-07-22.1 database, executed on an Amazon EC2 r6i.2xlarge instance with 8 parallel processes, yielding 349,758 new annotations.

All data were consolidated into a SQLite database (280 MB) containing nine tables: nuccore, taxonomy, typing, amr, plasmidfinder, typing_markers, pgap_features, biosample_location, and ncbi_extra. The database construction pipeline (build_db.py) is fully automated and downloads source files directly from Figshare, enabling reproducible rebuilds when PLSDB releases updates.

### Host source classification

Biosample entries were classified into six host/source categories (Human, Animal, Soil, Water, Food, Environment) using a hierarchical decision algorithm based on NCBI taxonomy IDs (9606 for *Homo sapiens*, livestock species for Animal), BioSample package types (e.g., "Pathogen: clinical or host-associated"), and ecosystem tags from the PLSDB biosample metadata. A total of 23,522 biosamples were classified, covering 57,751 geolocated plasmids.

### Machine learning analyses

**Mobility prediction**: An XGBoost classifier was trained on 57,580 plasmids with complete metadata (country, host genus, host source category, year, plasmid length, GC content, and binary indicators for the top 10 Inc groups) to predict mobility class (conjugative, mobilizable, non-mobilizable). Five-fold cross-validation was used to assess generalisation performance. SHAP TreeExplainer values were computed on a 5,000-plasmid subsample to quantify per-feature contributions to predictions (17).

**Co-mobilization prediction**: A Random Forest classifier was trained to predict whether a mobilizable plasmid shares a bacterial host (same BioSample accession) with a conjugative plasmid, using relaxase type, oriT type, Inc group, plasmid length, GC content, and host genus as features. Five-fold cross-validation assessed predictive accuracy.

**Integron carriage prediction**: A Random Forest classifier predicted the presence of class 1 integrons (defined as co-occurrence of qacEdelta1 and sul1) from plasmid features including mobility, genus, length, GC content, and Inc group indicators.

**UMAP embedding**: A two-dimensional UMAP projection (n_neighbors=30, min_dist=0.3) was computed on 15,000 randomly sampled plasmids using the encoded feature matrix, enabling visual identification of plasmid ecotypes (18).

### Confound decomposition

**Matched comparison**: Conjugative fractions were computed for each (species, country) pair with a minimum of 10 observations, controlling for species composition bias when comparing regional mobility patterns.

**Rarefaction analysis**: Each country was subsampled to equal sample sizes (50 to 2,000) with 50 bootstrap replicates per size level, computing the conjugative fraction and standard deviation at each level to assess sample size effects on observed regional differences.

**Simpson's paradox detection**: For each Inc group with more than 50 plasmids, the overall conjugative fraction was compared to per-species conjugative fractions. Cases with more than 20 percentage point divergence between species were flagged.

### Sequence analysis engine

A modular sequence analysis pipeline (seq_analysis.py) scans input DNA sequences for: (i) restriction enzyme recognition sites for 25 common Type II enzymes; (ii) sigma-70 promoter motifs (-35 TTGACA, spacer 14-20 bp, -10 TATAAT); (iii) Shine-Dalgarno ribosome binding site sequences; (iv) vector backbone signatures from common cloning vectors; (v) IS element terminal inverted repeats for 10 IS families; (vi) direct repeat/target site duplications (4-12 bp) flanking insertions; (vii) NTEKPC transposon structure around blaKPC genes; (viii) per-window codon adaptation index relative to *E. coli* K-12; and (ix) 4-mer frequency-based naturalness scoring. An engineering score (0-100) was computed using PLSDB-calibrated thresholds, with IS element presence reducing the score (natural mobile element markers) and synthetic-specific vector signatures (T7 promoter, lacZ alpha, f1 origin) increasing it. Vector backbone sequences shared between natural plasmids and synthetic constructs (ColE1, pBR322 origins) received reduced weighting.

### Retro-mobilization analysis

Non-mobilizable plasmids were assessed for potential horizontal transfer routes by identifying co-occurrence with conjugative plasmids in the same bacterial host (same BioSample UID), with mobilizable plasmids (relay mobilization), and with phage-compatible sizes (<10 kb for generalised transduction). Retro-mobilization — the ability of conjugative plasmids to mobilise DNA from recipient cells back into the donor — was quantified by counting non-mobilizable plasmids sharing a BioSample with at least one conjugative plasmid (19, 20).

### Web application and deployment

PlasmidNet was implemented using Dash 2.18.2 (21) with Plotly for interactive visualisation. The application is served via Gunicorn with a single worker process and four threads. BLAST+ is included via the heroku-community/apt buildpack for plasmid comparison features. The SQLite database is built during deployment by the bin/post_compile hook, which downloads source CSVs from Figshare and harvests additional plasmids from NCBI. A REST API implemented with Flask Blueprint provides programmatic access to all database queries. Source code is available at https://github.com/lcerdeira/plasmidnet under the MIT license.

## Results

### Database scope and integration

PlasmidNet integrates 182,362 complete bacterial plasmids: 72,556 from PLSDB 2025 with full MOB typing, AMR annotations, and taxonomic metadata, supplemented by 109,806 additional plasmids harvested from NCBI Nucleotide (Table 1). The combined AMR annotation set comprises 600,896 entries covering 163 drug classes. PGAP-derived annotations include 120,943 toxin-antitoxin system components across 40,759 plasmids, 54,343 phage-related elements on 15,340 plasmids, and 438,616 transposase/IS element annotations. Geographic metadata covers 57,751 plasmids from 100+ countries, with 23,522 classified by host source (Human: 11,464; Animal: 5,925; Soil: 1,576; Water: 1,332; Food: 1,059; Environment: 959).

### Plasmid mobility and transfer mechanisms

Of the 72,556 PLSDB plasmids with MOB typing, 22,588 (31.1%) were classified as conjugative, 21,152 (29.2%) as mobilizable, and 28,816 (39.7%) as non-mobilizable. Analysis of IS element distribution revealed IS26 as the most prevalent element (48,228 annotations across 13,298 plasmids), followed by IS3 (21,125 annotations) and IS5 (14,932 annotations). IS element family distribution varied significantly across Inc groups, with IncFIB/FII plasmids enriched for IS26 and IS1 elements.

Class 1 integrons (defined by co-occurrence of qacEdelta1 and sul1) were identified on 5,371 plasmids. Gene cassette analysis revealed sul1 (5,354), aadA2 (1,587), dfrA12 (1,227), arr-3 (983), catB3 (627), and blaOXA-1 (478) as the most frequent cassette-associated genes within 5 kb of qacEdelta1. Notably, 74% of integron-carrying plasmids were classified as conjugative, indicating high dissemination potential.

### Co-mobilization and retro-mobilization

Among mobilizable plasmids, 9,513 (45.0%) were found to co-occur in the same bacterial host as at least one conjugative plasmid, suggesting frequent co-mobilization via the conjugative T4SS. The co-mobilization ML predictor achieved 77.0% accuracy (5-fold CV), with relaxase type identified as the strongest predictor. Analysis of relaxase-T4SS compatibility revealed 885 co-occurrences of MOBC relaxase with MPF_F T4SS, a combination not described in standard compatibility rules — suggesting a potential novel mobilization pathway.

For non-mobilizable plasmids (n=28,816), multiple alternative transfer routes were identified: 9,555 (33.2%) shared a host with a conjugative plasmid (retro-mobilization potential), 10,736 (37.3%) had a mobilizable partner, and 8,585 (29.8%) were small enough (<10 kb) for phage-mediated transduction. Critically, 2,226 non-mobilizable plasmids carrying AMR genes co-existed with a conjugative partner, representing direct retro-mobilization risk for resistance dissemination.

### Machine learning for confound decomposition

The XGBoost mobility classifier achieved 82.6% accuracy on the full dataset. SHAP analysis revealed Inc group indicators (particularly IncFIB and IncFII) as the strongest predictors, followed by plasmid length, host genus, and GC content. Country of origin retained independent predictive power after controlling for these features, suggesting that regional differences in plasmid mobility are not fully explained by species composition or sequencing bias.

Rarefaction analysis confirmed that conjugative fractions stabilised at different values for different countries when subsampled to equal sizes, supporting the biological relevance of regional differences. No Simpson's paradoxes were detected — the direction of mobility differences was consistent across species within each country, though the magnitude varied.

### KPC transposon context

Analysis of 4,800+ blaKPC-carrying plasmids revealed blaTEM as the most common co-located gene within 5 kb (192 plasmids), consistent with the NTEKPC-IId element structure. The majority (59%) of KPC plasmids were conjugative, predominantly on IncFIB/FII (204 plasmids), IncFIA/FII/IncR (168), and IncN (161) backbones. The conserved genetic context of blaKPC across diverse plasmid backbones and host species — including flanking IS elements, target site duplications, and consistent gene arrangements — supports the classification of NTEKPC elements as bona fide transposons rather than arbitrary genetic islands.

### Geographic and temporal patterns

Temporal analysis (2010-2024) revealed increasing prevalence of several resistance classes, with beta-lactam and aminoglycoside resistance showing the steepest rises. Inc group trends showed expansion of IncFIB and IncFII plasmids relative to other replicon types. Geographic comparison demonstrated significant variation in both AMR profiles and mobility patterns across countries, with host source (human vs. animal vs. environmental) explaining a substantial portion of the variance.

### Sequence analysis and engineering detection

The sequence analysis engine was validated against known natural and engineered plasmids. The natural KPC-carrying plasmid MH595533 scored 45/100 (Ambiguous), correctly reflecting its use of ColE1-type replication (shared between natural IncQ plasmids and laboratory vectors) while lacking synthetic-specific signatures. The engineered vector pUC19 scored 67/100 (Likely engineered), correctly identified by its f1 origin and lacZ alpha — signatures unique to synthetic constructs. The IS element-aware scoring reduced false positive engineering classifications for natural resistance platforms that accumulate restriction site clusters through transposon insertion rather than deliberate cloning.

## Discussion

PlasmidNet addresses a critical gap in the plasmid genomics toolkit by providing an integrated platform that combines the largest unified plasmid dataset with interactive analytics, machine learning, and programmatic access. Several aspects of the platform warrant discussion.

The identification of retro-mobilization potential for one-third of non-mobilizable plasmids has important implications for AMR surveillance. While retro-mobilization has been demonstrated experimentally in *Bacillus* systems (20), its prevalence in clinical Enterobacteriaceae has not been systematically assessed. Our finding that 2,226 AMR-carrying non-mobilizable plasmids share hosts with conjugative plasmids suggests that resistance gene mobilisation may occur more frequently than predicted by mobility classification alone.

The co-mobilization analysis revealed an unexpected MOBC-MPF_F compatibility, representing 885 co-occurrences not described in standard relaxase-T4SS compatibility rules (7). This finding warrants experimental validation, as it could represent a previously unrecognised mobilization pathway or reflect indirect co-selection rather than functional T4SS interaction.

The engineering detection module addresses the growing need for forensic genomics tools that can distinguish natural from synthetic plasmid sequences (22). By calibrating against the PLSDB natural plasmid baseline and incorporating IS element detection, PlasmidNet reduces false positive classifications of natural resistance platforms — a limitation of approaches that rely solely on restriction site density or codon usage bias.

Several limitations should be noted. The NCBI-harvested plasmids (109,806 entries) currently lack MOB typing and PGAP annotations, limiting their use in mobility and IS element analyses. The batch-mode MOBsuite typing used in the EC2 pipeline produced per-batch rather than per-sequence results for these entries. Future versions will incorporate individual plasmid typing. The sequence analysis engineering score, while validated against exemplar sequences, has not been systematically benchmarked against a large dataset of confirmed natural and engineered plasmids. The geographic analyses are subject to sequencing bias, as different countries have different surveillance priorities and sequencing capacity — a limitation partially addressed by our rarefaction and confound decomposition approaches.

In conclusion, PlasmidNet provides a comprehensive, freely accessible platform for plasmid analytics that integrates database exploration, machine learning, and sequence analysis. By combining 182,362 plasmids with interactive visualisation and a REST API, it enables both hypothesis-driven research and exploratory analysis of plasmid-mediated AMR dissemination.

## Data availability

The PlasmidNet database, source code, and documentation are freely available:
- Dashboard: https://plasmidnet-7302570f4818.herokuapp.com/
- Source code: https://github.com/lcerdeira/plasmidnet (MIT license)
- Documentation: https://plasmidnet.readthedocs.io/
- REST API: https://plasmidnet-7302570f4818.herokuapp.com/api/

## Acknowledgements

The author acknowledges the PLSDB team at Saarland University for maintaining the plasmid database and providing open access to their data through Figshare.

## Funding

This work was supported by [funding information to be added].

## Conflict of Interest

None declared.

## References

1. Partridge SR, Kwong SM, Firth N, Jensen SO. Mobile genetic elements associated with antimicrobial resistance. *Clin Microbiol Rev*. 2018;31(4):e00088-17.

2. Carattoli A. Plasmids and the spread of resistance. *Int J Med Microbiol*. 2013;303(6-7):298-304.

3. Tacconelli E, Carrara E, Savoldi A, et al. Discovery, research, and development of new antibiotics: the WHO priority list of antibiotic-resistant bacteria and tuberculosis. *Lancet Infect Dis*. 2018;18(3):318-327.

4. Schmartz GP, Schmitt K, Backes C, et al. PLSDB: advancing a comprehensive database of bacterial plasmids. *Nucleic Acids Res*. 2025;53(D1):D235-D241.

5. Galata V, Fehlmann T, Backes C, Keller A. PLSDB: a resource of complete bacterial plasmids. *Nucleic Acids Res*. 2019;47(D1):D195-D202.

6. Carattoli A, Zankari E, Garcia-Fernandez A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother*. 2014;58(7):3895-3903.

7. Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microb Genom*. 2018;4(8):e000206.

8. Alcock BP, Raphenya AR, Lau TTY, et al. CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. *Nucleic Acids Res*. 2020;48(D1):D482-D489.

9. Feldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep*. 2021;11(1):12728.

10. Siguier P, Gourbeyre E, Chandler M. Bacterial insertion sequences: their genomic impact and diversity. *FEMS Microbiol Rev*. 2014;38(5):865-891.

11. He S, Hickman AB, Varani AM, et al. Insertion sequence IS26 reorganizes plasmids in clinically isolated multidrug-resistant bacteria by replicative transposition. *mBio*. 2015;6(3):e00762.

12. Gillings MR. Integrons: past, present, and future. *Microbiol Mol Biol Rev*. 2014;78(2):257-277.

13. Naas T, Cuzon G, Villegas MV, Lartigue MF, Quinn JP, Nordmann P. Genetic structures at the origin of acquisition of the beta-lactamase blaKPC gene. *Antimicrob Agents Chemother*. 2008;52(4):1257-1263.

14. Chen L, Chavda KD, Melano RG, et al. Comparative genomic analysis of KPC-encoding pKpQIL-like plasmids and their distribution in New Jersey and New York hospitals. *Antimicrob Agents Chemother*. 2014;58(5):2871-2877.

15. Arango-Argoty G, Garner E, Pruden A, Heath LS, Vikesland P, Zhang L. DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. *Microbiome*. 2018;6(1):23.

16. Ruiz-Perez D, Guan H, Madhivanan P, et al. So you think you can PLS-DA? *BMC Bioinformatics*. 2020;21(1):2.

17. Lundberg SM, Erion G, Chen H, et al. From local explanations to global understanding with explainable AI for trees. *Nat Mach Intell*. 2020;2(1):56-67.

18. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv*. 2018;1802.03426.

19. Ramsay JP, Firth N. Diverse mobilization strategies facilitate transfer of non-conjugative mobile genetic elements. *Curr Opin Microbiol*. 2017;38:1-9.

20. Timmery S, Hu X, Mahillon J. Characterization of bacilli isolated from the confined environments of the Antarctic Concordia Station and the International Space Station. *Astrobiology*. 2011;11(4):323-334.

21. Plotly Technologies Inc. Dash: A productive Python framework for building web analytic applications. 2024. https://dash.plotly.com/

22. Dunlap G, Pauwels E. The Toolbox: biosecurity and AI. *Johns Hopkins Center for Health Security*. 2023.

23. Douarre PE, Mallet L, Raber N, Felten A, Feurer C. Analysis of COMPASS, a new comprehensive plasmid database revealed prevalence of multireplicon and extensive diversity of IncF plasmids. *Front Microbiol*. 2020;11:483.

24. Li X, Xie Y, Liu M, et al. oriTfinder: a web-based tool for the identification of origin of transfers in DNA sequences of bacterial mobile genetic elements. *Nucleic Acids Res*. 2018;46(W1):W229-W234.

25. Redondo-Salvo S, Fernandez-Lopez R, Ruiz R, et al. Pathways for horizontal gene transfer in bacteria revealed by a global map of their plasmids. *Nat Commun*. 2020;11(1):3602.

26. Acman M, van Dorp L, Balloux F. Large-scale network analysis captures biological features of bacterial plasmids. *Nat Commun*. 2020;11(1):2452.

27. Shintani M, Sanchez ZK, Kimbara K. Genomics of microbial plasmids: classification and identification based on replication and transfer systems and host taxonomy. *Front Microbiol*. 2015;6:242.
