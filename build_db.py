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
            PMLST_alleles     TEXT
        )
    """)
    with open(os.path.join(RAW_DIR, "typing.csv"), newline="") as f:
        reader = csv.DictReader(f)
        rows = 0
        for row in reader:
            try:
                cur.execute("""
                    INSERT OR IGNORE INTO typing VALUES (?,?,?,?,?,?,?,?,?,?,?)
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
    conn.commit()

    # ── Stats ───────────────────────────────────────────────────────
    print("\n=== Database Summary ===")
    for table in ["nuccore", "taxonomy", "typing", "amr", "plasmidfinder", "typing_markers"]:
        count = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table:20s} {count:>10,} rows")

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
