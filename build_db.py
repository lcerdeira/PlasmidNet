#!/usr/bin/env python3
"""
Build a SQLite database from PLSDB Figshare CSV/TSV files.

Usage:
    python build_db.py

Reads from data/raw/*.csv and produces data/plsdb.db (~50 MB).
"""

import csv
import os
import sqlite3
import sys

RAW_DIR = os.path.join(os.path.dirname(__file__), "data", "raw")
DB_PATH = os.path.join(os.path.dirname(__file__), "data", "plsdb.db")


def main():
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        print(f"Removed old {DB_PATH}")

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    cur = conn.cursor()

    # ── 1. nuccore (plasmid metadata) ───────────────────────────────
    print("Loading nuccore.csv ...")
    cur.execute("""
        CREATE TABLE nuccore (
            NUCCORE_UID       INTEGER,
            NUCCORE_ACC       TEXT PRIMARY KEY,
            NUCCORE_Description TEXT,
            NUCCORE_CreateDate TEXT,
            NUCCORE_Completeness TEXT,
            NUCCORE_Length    INTEGER,
            NUCCORE_Source    TEXT,
            BIOSAMPLE_UID    INTEGER,
            TAXONOMY_UID     INTEGER,
            NUCCORE_GC       REAL,
            NUCCORE_Topology  TEXT
        )
    """)
    with open(os.path.join(RAW_DIR, "nuccore.csv"), newline="") as f:
        reader = csv.DictReader(f)
        rows = 0
        for row in reader:
            try:
                cur.execute("""
                    INSERT OR IGNORE INTO nuccore VALUES (?,?,?,?,?,?,?,?,?,?,?)
                """, (
                    _int(row.get("NUCCORE_UID")),
                    row.get("NUCCORE_ACC", ""),
                    row.get("NUCCORE_Description", ""),
                    row.get("NUCCORE_CreateDate", ""),
                    row.get("NUCCORE_Completeness", ""),
                    _int(row.get("NUCCORE_Length")),
                    row.get("NUCCORE_Source", ""),
                    _int(row.get("BIOSAMPLE_UID")),
                    _int(row.get("TAXONOMY_UID")),
                    _float(row.get("NUCCORE_GC")),
                    row.get("NUCCORE_Topology", ""),
                ))
                rows += 1
            except Exception as e:
                pass
    conn.commit()
    print(f"  {rows:,} plasmids loaded")

    # ── 2. taxonomy ─────────────────────────────────────────────────
    print("Loading taxonomy.csv ...")
    cur.execute("""
        CREATE TABLE taxonomy (
            TAXONOMY_UID     INTEGER PRIMARY KEY,
            TAXONOMY_superkingdom TEXT,
            TAXONOMY_phylum  TEXT,
            TAXONOMY_class   TEXT,
            TAXONOMY_order   TEXT,
            TAXONOMY_family  TEXT,
            TAXONOMY_genus   TEXT,
            TAXONOMY_species TEXT,
            TAXONOMY_strain  TEXT
        )
    """)
    with open(os.path.join(RAW_DIR, "taxonomy.csv"), newline="") as f:
        reader = csv.DictReader(f)
        rows = 0
        for row in reader:
            try:
                cur.execute("""
                    INSERT OR IGNORE INTO taxonomy VALUES (?,?,?,?,?,?,?,?,?)
                """, (
                    _int(row.get("TAXONOMY_UID")),
                    row.get("TAXONOMY_superkingdom", ""),
                    row.get("TAXONOMY_phylum", ""),
                    row.get("TAXONOMY_class", ""),
                    row.get("TAXONOMY_order", ""),
                    row.get("TAXONOMY_family", ""),
                    row.get("TAXONOMY_genus", ""),
                    row.get("TAXONOMY_species", ""),
                    row.get("TAXONOMY_strain", ""),
                ))
                rows += 1
            except Exception:
                pass
    conn.commit()
    print(f"  {rows:,} taxonomy entries")

    # ── 3. typing (MOB typing, Inc groups, mobility, pMLST) ────────
    print("Loading typing.csv ...")
    cur.execute("""
        CREATE TABLE typing (
            NUCCORE_ACC       TEXT PRIMARY KEY,
            rep_type          TEXT,
            relaxase_type     TEXT,
            mpf_type          TEXT,
            orit_type         TEXT,
            predicted_mobility TEXT,
            predicted_host_range_overall_name TEXT,
            predicted_host_range_overall_rank TEXT,
            PMLST_scheme      TEXT,
            PMLST_sequence_type TEXT,
            PMLST_alleles     TEXT,
            associated_pmid   TEXT
        )
    """)
    with open(os.path.join(RAW_DIR, "typing.csv"), newline="") as f:
        reader = csv.DictReader(f)
        rows = 0
        for row in reader:
            try:
                cur.execute("""
                    INSERT OR IGNORE INTO typing VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
                """, (
                    row.get("NUCCORE_ACC", ""),
                    row.get("rep_type(s)", ""),
                    row.get("relaxase_type(s)", ""),
                    row.get("mpf_type", ""),
                    row.get("orit_type(s)", ""),
                    row.get("predicted_mobility", ""),
                    row.get("predicted_host_range_overall_name", ""),
                    row.get("predicted_host_range_overall_rank", ""),
                    row.get("PMLST_scheme", ""),
                    row.get("PMLST_sequence_type", ""),
                    row.get("PMLST_alleles", ""),
                    row.get("associated_pmid(s)", ""),
                ))
                rows += 1
            except Exception:
                pass
    conn.commit()
    print(f"  {rows:,} typing entries")

    # ── 4. amr (antimicrobial resistance genes) ─────────────────────
    print("Loading amr.tsv ...")
    cur.execute("""
        CREATE TABLE amr (
            NUCCORE_ACC       TEXT,
            gene_symbol       TEXT,
            gene_name         TEXT,
            drug_class        TEXT,
            antimicrobial_agent TEXT,
            input_gene_start  INTEGER,
            input_gene_stop   INTEGER,
            strand_orientation TEXT,
            sequence_identity REAL,
            coverage_percentage REAL
        )
    """)
    with open(os.path.join(RAW_DIR, "amr.tsv"), newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = 0
        batch = []
        for row in reader:
            batch.append((
                row.get("NUCCORE_ACC", ""),
                row.get("gene_symbol", ""),
                row.get("gene_name", ""),
                row.get("drug_class", ""),
                row.get("antimicrobial_agent", ""),
                _int(row.get("input_gene_start")),
                _int(row.get("input_gene_stop")),
                row.get("strand_orientation", ""),
                _float(row.get("sequence_identity")),
                _float(row.get("coverage_percentage")),
            ))
            rows += 1
            if len(batch) >= 5000:
                cur.executemany("INSERT INTO amr VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
                batch = []
        if batch:
            cur.executemany("INSERT INTO amr VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
    conn.commit()
    print(f"  {rows:,} AMR annotations")

    # ── 5. plasmidfinder ────────────────────────────────────────────
    print("Loading plasmidfinder.csv ...")
    cur.execute("""
        CREATE TABLE plasmidfinder (
            NUCCORE_ACC TEXT,
            typing      TEXT,
            identity    REAL,
            coverage    REAL,
            organism_L1 TEXT,
            organism_L2 TEXT
        )
    """)
    with open(os.path.join(RAW_DIR, "plasmidfinder.csv"), newline="") as f:
        reader = csv.DictReader(f)
        rows = 0
        batch = []
        for row in reader:
            batch.append((
                row.get("NUCCORE_ACC", ""),
                row.get("typing", ""),
                _float(row.get("identity")),
                _float(row.get("coverage")),
                row.get("organism_L1", ""),
                row.get("organism_L2", ""),
            ))
            rows += 1
            if len(batch) >= 5000:
                cur.executemany("INSERT INTO plasmidfinder VALUES (?,?,?,?,?,?)", batch)
                batch = []
        if batch:
            cur.executemany("INSERT INTO plasmidfinder VALUES (?,?,?,?,?,?)", batch)
    conn.commit()
    print(f"  {rows:,} plasmidfinder entries")

    # ── 6. typing_markers (MOB detail markers with positions) ──────
    print("Loading typing_markers.csv ...")
    cur.execute("""
        CREATE TABLE typing_markers (
            NUCCORE_ACC TEXT,
            element     TEXT,
            sstart      INTEGER,
            send        INTEGER,
            sstrand     TEXT,
            pident      REAL,
            qcovhsp     REAL,
            biomarker   TEXT
        )
    """)
    with open(os.path.join(RAW_DIR, "typing_markers.csv"), newline="") as f:
        reader = csv.DictReader(f)
        rows = 0
        batch = []
        for row in reader:
            batch.append((
                row.get("NUCCORE_ACC", ""),
                row.get("element", ""),
                _int(row.get("sstart")),
                _int(row.get("send")),
                row.get("sstrand", ""),
                _float(row.get("pident")),
                _float(row.get("qcovhsp")),
                row.get("biomarker", row.get("element", "")),
            ))
            rows += 1
            if len(batch) >= 5000:
                cur.executemany("INSERT INTO typing_markers VALUES (?,?,?,?,?,?,?,?)", batch)
                batch = []
        if batch:
            cur.executemany("INSERT INTO typing_markers VALUES (?,?,?,?,?,?,?,?)", batch)
    conn.commit()
    print(f"  {rows:,} typing marker entries")

    # ── 7. biosample_location (geographic data) ──────────────────────
    bio_path = os.path.join(RAW_DIR, "biosample.csv")
    if os.path.exists(bio_path):
        print("Loading biosample.csv (locations) ...")
        cur.execute("""
            CREATE TABLE biosample_location (
                BIOSAMPLE_UID INTEGER PRIMARY KEY,
                lat REAL, lng REAL,
                location_raw TEXT, country TEXT, ecosystem TEXT,
                host_category TEXT
            )
        """)

        COUNTRY_BOUNDS = {
            "USA": (24, 50, -125, -66), "China": (18, 54, 73, 135),
            "Japan": (24, 46, 123, 146), "South Korea": (33, 39, 124, 132),
            "India": (6, 36, 68, 98), "United Kingdom": (49, 61, -11, 2),
            "Germany": (47, 55, 5, 16), "France": (41, 51, -5, 10),
            "Spain": (36, 44, -10, 5), "Brazil": (-34, 6, -74, -34),
            "Australia": (-44, -10, 113, 154), "Canada": (42, 72, -141, -52),
            "Taiwan": (21, 26, 119, 123), "Thailand": (5, 21, 97, 106),
            "Myanmar": (9, 29, 92, 102), "Bangladesh": (20, 27, 88, 93),
            "Denmark": (54, 58, 8, 16), "Switzerland": (45, 48, 5, 11),
            "Netherlands": (50, 54, 3, 8), "Italy": (36, 47, 6, 19),
            "Sweden": (55, 70, 10, 25), "Russia": (41, 82, 19, 180),
            "Mexico": (14, 33, -118, -86), "Vietnam": (8, 24, 102, 110),
            "Indonesia": (-11, 6, 95, 141), "Singapore": (1, 2, 103, 104),
            "Israel": (29, 34, 34, 36), "Turkey": (35, 42, 25, 45),
            "Colombia": (-5, 14, -82, -66), "Nigeria": (4, 14, 2, 15),
            "South Africa": (-35, -22, 16, 33), "Kenya": (-5, 6, 33, 42),
            "Egypt": (22, 32, 24, 37), "Pakistan": (23, 37, 60, 78),
            "Iran": (25, 40, 44, 64), "Argentina": (-55, -21, -74, -53),
            "Norway": (57, 72, 4, 32), "Poland": (49, 55, 14, 25),
            "Belgium": (49, 52, 2, 7), "Portugal": (36, 42, -10, -6),
            "Greece": (34, 42, 19, 30), "Saudi Arabia": (16, 33, 34, 56),
            "Philippines": (4, 21, 116, 127), "Malaysia": (0, 8, 99, 120),
            "New Zealand": (-47, -34, 166, 179),
        }
        COUNTRY_FIXES = {
            "United States": "USA", "United States of America": "USA",
            "Republic of Korea": "South Korea", "Korea": "South Korea",
            "United Kingdom of Great Britain and Northern Ireland": "United Kingdom",
            "Viet Nam": "Vietnam", "Russian Federation": "Russia",
            "Brasil": "Brazil", "Czech Republic": "Czechia",
        }

        def _normalize_country(raw, lat, lng):
            if not raw:
                return ""
            parts = raw.split(",")
            country = parts[0].strip()
            if len(country) > 25:
                # Try lat/lng lookup
                for c, (la1, la2, lo1, lo2) in COUNTRY_BOUNDS.items():
                    if la1 <= lat <= la2 and lo1 <= lng <= lo2:
                        return c
                return ""
            return COUNTRY_FIXES.get(country, country)

        HUMAN_TAXIDS = {9606}
        ANIMAL_TAXIDS = {9823, 9913, 9031, 9615, 9685, 9940, 9796, 8839,
                         9825, 10090, 7460, 8030, 9986, 9598, 9544}

        def _classify_host(eco, pkg, taxid):
            eco_l = eco.lower()
            pkg_l = pkg.lower()
            if taxid in HUMAN_TAXIDS:
                return "Human"
            if taxid in ANIMAL_TAXIDS:
                return "Animal"
            if "clinical" in pkg_l or "human-associated" in pkg_l:
                return "Human"
            if any(k in eco_l for k in ("blood", "urinary", "respiratory",
                                         "circulatory")):
                return "Human"
            if any(k in eco_l for k in ("fecal", "gastrointestinal", "rectal")):
                return "Human" if taxid in HUMAN_TAXIDS or taxid < 0 else "Animal"
            if "host_associated" in eco_l:
                return "Human" if taxid in HUMAN_TAXIDS or taxid < 0 else "Animal"
            if any(k in eco_l for k in ("soil", "terrestrial", "rhizosphere")):
                return "Soil"
            if any(k in eco_l for k in ("aquatic", "water", "wastewater",
                                         "marine", "freshwater")):
                return "Water"
            if any(k in eco_l for k in ("food", "meat", "dairy", "fermented",
                                         "drink", "milk")):
                return "Food"
            if any(k in eco_l for k in ("farm", "manure", "livestock")):
                return "Animal"
            if any(k in eco_l for k in ("hospital", "anthropogenic")):
                return "Environment"
            if any(k in pkg_l for k in ("environmental", "wastewater")):
                return "Environment"
            if "water" in pkg_l:
                return "Water"
            if "soil" in pkg_l:
                return "Soil"
            return "Other"

        batch = []
        rows = 0
        with open(bio_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                buid = row.get("BIOSAMPLE_UID", "")
                try:
                    buid = int(float(buid))
                except (ValueError, TypeError):
                    continue
                lat_s = row.get("LOCATION_lat", "").strip()
                lng_s = row.get("LOCATION_lng", "").strip()
                loc_raw = row.get("LOCATION_name", "").strip()
                eco = row.get("ECOSYSTEM_tags", "").strip()
                pkg = row.get("BIOSAMPLE_package", "").strip()
                taxid = -1
                try:
                    taxid = int(float(row.get("ECOSYSTEM_taxid", -1)))
                except (ValueError, TypeError):
                    pass
                if lat_s and lng_s and loc_raw:
                    try:
                        lat_f, lng_f = float(lat_s), float(lng_s)
                    except ValueError:
                        continue
                    country = _normalize_country(loc_raw, lat_f, lng_f)
                    host_cat = _classify_host(eco, pkg, taxid)
                    if country:
                        batch.append((buid, lat_f, lng_f, loc_raw, country,
                                      eco, host_cat))
                        rows += 1
                        if len(batch) >= 5000:
                            cur.executemany(
                                "INSERT OR IGNORE INTO biosample_location "
                                "VALUES (?,?,?,?,?,?,?)", batch)
                            batch = []
        if batch:
            cur.executemany(
                "INSERT OR IGNORE INTO biosample_location VALUES (?,?,?,?,?,?,?)",
                batch)
        cur.execute("CREATE INDEX idx_bsloc_country ON biosample_location(country)")
        cur.execute("CREATE INDEX idx_bsloc_host ON biosample_location(host_category)")
        conn.commit()
        print(f"  {rows:,} biosample locations")
    else:
        print("biosample.csv not found, skipping location data")

    # ── 8. pgap_features (TA, phage, transposase from proteins.csv) ─
    pgap_path = os.path.join(RAW_DIR, "pgap_filtered.csv")
    if os.path.exists(pgap_path):
        print("Loading pgap_filtered.csv ...")
        cur.execute("""
            CREATE TABLE pgap_features (
                NUCCORE_ACC TEXT,
                gene TEXT,
                product TEXT,
                locus_tag TEXT,
                category TEXT,
                start INTEGER,
                end INTEGER,
                strand INTEGER
            )
        """)
        TA_GENES = {"rele", "relb", "higa", "higb", "mazf", "maze", "ccda", "ccdb",
                    "pard", "pare", "vapb", "vapc", "phd", "doc", "yoeb", "yefm",
                    "mqsr", "mqsa", "hipa", "hipb", "pema", "pemk", "pemi",
                    "hicb", "hica", "yafq", "dina", "brna", "brnb"}
        TA_KW = ["toxin-antitoxin", "antitoxin", "addiction module",
                 "plasmid stabilization", "killer protein",
                 "post-segregational", "type ii toxin"]
        PHAGE_KW = ["phage", "prophage", "bacteriophage", "tail fiber",
                    "capsid", "terminase", "portal protein", "baseplate",
                    "holin", "lysozyme"]
        with open(pgap_path, newline="") as f:
            reader = csv.DictReader(f)
            rows = 0
            batch = []
            for row in reader:
                gene = (row.get("gene") or "").lower()
                product = (row.get("product") or "").lower()
                if gene in TA_GENES or any(kw in product for kw in TA_KW):
                    cat = "TA"
                elif any(kw in product for kw in PHAGE_KW):
                    cat = "PHAGE"
                else:
                    cat = "TRANSPOSASE"
                strand = 1
                try:
                    strand = int(float(row.get("strand", 1)))
                except (ValueError, TypeError):
                    pass
                batch.append((
                    row.get("NUCCORE_ACC", ""),
                    row.get("gene") or row.get("locus_tag") or "",
                    row.get("product", ""),
                    row.get("locus_tag", ""),
                    cat,
                    _int(row.get("start")),
                    _int(row.get("end")),
                    strand,
                ))
                rows += 1
                if len(batch) >= 10000:
                    cur.executemany("INSERT INTO pgap_features VALUES (?,?,?,?,?,?,?,?)", batch)
                    batch = []
            if batch:
                cur.executemany("INSERT INTO pgap_features VALUES (?,?,?,?,?,?,?,?)", batch)
        conn.commit()
        print(f"  {rows:,} PGAP feature entries")
    else:
        print("pgap_filtered.csv not found, skipping PGAP features")

    # ── Indexes ─────────────────────────────────────────────────────
    print("Creating indexes ...")
    cur.execute("CREATE INDEX idx_amr_acc ON amr(NUCCORE_ACC)")
    cur.execute("CREATE INDEX idx_amr_drug ON amr(drug_class)")
    cur.execute("CREATE INDEX idx_amr_gene ON amr(gene_symbol)")
    cur.execute("CREATE INDEX idx_typing_mobility ON typing(predicted_mobility)")
    cur.execute("CREATE INDEX idx_typing_rep ON typing(rep_type)")
    cur.execute("CREATE INDEX idx_nuccore_topo ON nuccore(NUCCORE_Topology)")
    cur.execute("CREATE INDEX idx_nuccore_tax ON nuccore(TAXONOMY_UID)")
    cur.execute("CREATE INDEX idx_pf_acc ON plasmidfinder(NUCCORE_ACC)")
    cur.execute("CREATE INDEX idx_tm_acc ON typing_markers(NUCCORE_ACC)")
    if os.path.exists(pgap_path):
        cur.execute("CREATE INDEX idx_pgap_acc ON pgap_features(NUCCORE_ACC)")
        cur.execute("CREATE INDEX idx_pgap_cat ON pgap_features(category)")
    conn.commit()

    # ── 9. ncbi_extra (plasmids from NCBI not in PLSDB) ──────────────
    ncbi_path = os.path.join(RAW_DIR, "ncbi_plasmids.csv")
    if os.path.exists(ncbi_path):
        print("Loading ncbi_plasmids.csv ...")
        cur.execute("""
            CREATE TABLE ncbi_extra (
                accession    TEXT PRIMARY KEY,
                description  TEXT,
                length       INTEGER,
                topology     TEXT,
                create_date  TEXT,
                organism     TEXT,
                taxid        INTEGER
            )
        """)
        with open(ncbi_path, newline="") as f:
            reader = csv.DictReader(f)
            rows = 0
            batch = []
            for row in reader:
                batch.append((
                    row.get("accession", ""),
                    row.get("description", ""),
                    _int(row.get("length")),
                    row.get("topology", ""),
                    row.get("create_date", ""),
                    row.get("organism", ""),
                    _int(row.get("taxid")),
                ))
                rows += 1
                if len(batch) >= 5000:
                    cur.executemany(
                        "INSERT OR IGNORE INTO ncbi_extra VALUES (?,?,?,?,?,?,?)",
                        batch)
                    batch = []
            if batch:
                cur.executemany(
                    "INSERT OR IGNORE INTO ncbi_extra VALUES (?,?,?,?,?,?,?)",
                    batch)
        cur.execute("CREATE INDEX idx_ncbi_date ON ncbi_extra(create_date)")
        cur.execute("CREATE INDEX idx_ncbi_org ON ncbi_extra(organism)")
        conn.commit()
        print(f"  {rows:,} NCBI extra plasmids")
    else:
        print("ncbi_plasmids.csv not found - run ncbi_harvest.py first")

    # ── Stats ───────────────────────────────────────────────────────
    print("\n=== Database Summary ===")
    tables = ["nuccore", "taxonomy", "typing", "amr", "plasmidfinder",
              "typing_markers"]
    if os.path.exists(pgap_path):
        tables.append("pgap_features")
    tables.append("biosample_location")
    if os.path.exists(ncbi_path):
        tables.append("ncbi_extra")
    for table in tables:
        try:
            count = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            print(f"  {table:20s} {count:>10,} rows")
        except Exception:
            pass

    # Total unique plasmids
    plsdb_count = cur.execute("SELECT COUNT(*) FROM nuccore").fetchone()[0]
    ncbi_count = 0
    try:
        ncbi_count = cur.execute("SELECT COUNT(*) FROM ncbi_extra").fetchone()[0]
    except Exception:
        pass
    print(f"\n  Total unique plasmids: {plsdb_count + ncbi_count:,} "
          f"(PLSDB {plsdb_count:,} + NCBI {ncbi_count:,})")

    conn.close()
    size_mb = os.path.getsize(DB_PATH) / 1e6
    print(f"\nDatabase: {DB_PATH} ({size_mb:.1f} MB)")
    print("Done!")


def _int(val):
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return None


def _float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    main()
