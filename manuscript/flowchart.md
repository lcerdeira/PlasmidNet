# PlasmidNet Architecture Flowchart

```
┌─────────────────────────────────────────────────────────────────────┐
│                        DATA SOURCES                                 │
├──────────────────────┬──────────────────────┬───────────────────────┤
│   PLSDB 2025         │   NCBI Nucleotide    │   NCBI GenBank        │
│   (Figshare)         │   (E-utilities)      │   (on-demand)         │
│                      │                      │                       │
│  • nuccore.csv       │  • ncbi_harvest.py   │  • efetch per         │
│  • typing.csv        │  • 109,806 complete  │    accession          │
│  • amr.tsv           │    plasmids          │  • CDS features       │
│  • taxonomy.csv      │  • Metadata only     │  • DNA sequence       │
│  • biosample.csv     │                      │  • GC content         │
│  • proteins.csv      │                      │                       │
│    (filtered→PGAP)   │                      │                       │
│  • typing_markers    │                      │                       │
│  • plasmidfinder     │                      │                       │
├──────────────────────┴──────────────────────┴───────────────────────┤
│                                                                     │
│                    ┌──────────────────────┐                          │
│                    │    build_db.py       │                          │
│                    │  SQLite Database     │                          │
│                    │    (280 MB)          │                          │
│                    └────────┬─────────────┘                          │
│                             │                                        │
├─────────────────────────────┼────────────────────────────────────────┤
│                     DATABASE LAYER                                   │
│                             │                                        │
│  ┌──────────────────────────┼──────────────────────────────────┐     │
│  │              db.py (SQLite queries)                         │     │
│  │                                                            │     │
│  │  ┌─────────┐ ┌─────────┐ ┌──────────┐ ┌───────────────┐  │     │
│  │  │ nuccore │ │ typing  │ │   amr    │ │ pgap_features │  │     │
│  │  │ 72,556  │ │ 72,556  │ │ 600,896  │ │   613,902     │  │     │
│  │  └─────────┘ └─────────┘ └──────────┘ └───────────────┘  │     │
│  │  ┌─────────┐ ┌──────────────┐ ┌──────────┐ ┌──────────┐  │     │
│  │  │taxonomy │ │biosample_loc │ │ncbi_extra│ │plasmidfndr│  │     │
│  │  │  7,439  │ │   23,522     │ │ 109,806  │ │  60,619   │  │     │
│  │  └─────────┘ └──────────────┘ └──────────┘ └──────────┘  │     │
│  └────────────────────────────────────────────────────────────┘     │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│                     ANALYSIS ENGINES                                 │
│                                                                     │
│  ┌────────────────┐  ┌─────────────────┐  ┌──────────────────┐     │
│  │ seq_analysis.py│  │plasmid_compare.py│  │   api.py         │     │
│  │                │  │                 │  │                  │     │
│  │ • RE sites     │  │ • BLASTn align  │  │ • REST endpoints │     │
│  │ • Promoters    │  │ • pyCirclize    │  │ • JSON responses │     │
│  │ • RBS          │  │ • GBK/FA export │  │ • Search/filter  │     │
│  │ • IS elements  │  │                 │  │                  │     │
│  │ • Direct rpts  │  └─────────────────┘  └──────────────────┘     │
│  │ • NTEKPC       │                                                 │
│  │ • Codon bias   │  ┌─────────────────────────────────────────┐   │
│  │ • K-mer score  │  │        ML Models (scikit-learn,         │   │
│  │ • Eng. score   │  │         XGBoost, SHAP, UMAP)            │   │
│  │ • Retro-mob    │  │                                         │   │
│  └────────────────┘  │  • Mobility predictor (82.6% acc)      │   │
│                      │  • Co-mobilization predictor (77% acc)  │   │
│                      │  • Integron carriage predictor           │   │
│                      │  • UMAP clustering (15K plasmids)       │   │
│                      │  • Rarefaction bootstrap                │   │
│                      │  • Simpson's paradox detection          │   │
│                      └─────────────────────────────────────────┘   │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│                     WEB APPLICATION                                  │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────────┐    │
│  │                 app.py (Dash + Plotly)                       │    │
│  │                                                             │    │
│  │  ┌─────────┐ ┌────────┐ ┌─────┐ ┌────────┐ ┌───────────┐  │    │
│  │  │Overview │ │Taxonomy│ │ AMR │ │ Viewer │ │Inc/MOB/IS │  │    │
│  │  └─────────┘ └────────┘ └─────┘ └────────┘ └───────────┘  │    │
│  │  ┌────────────┐ ┌─────────┐ ┌─────────┐ ┌───────┐        │    │
│  │  │Correlations│ │Geography│ │Analytics│ │Compare│        │    │
│  │  └────────────┘ └─────────┘ └─────────┘ └───────┘        │    │
│  │  ┌────────────┐ ┌─────┐                                   │    │
│  │  │Seq Analysis│ │ API │                                   │    │
│  │  └────────────┘ └─────┘                                   │    │
│  └─────────────────────────────────────────────────────────────┘    │
│                             │                                        │
│                    ┌────────┴────────┐                               │
│                    │   Deployment    │                               │
│                    │                 │                               │
│                    │  • Heroku       │                               │
│                    │  • Gunicorn     │                               │
│                    │  • BLAST+ (Apt)│                               │
│                    └─────────────────┘                               │
└─────────────────────────────────────────────────────────────────────┘
```
