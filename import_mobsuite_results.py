#!/usr/bin/env python3
"""
Import MOBsuite + AMRFinderPlus results from EC2 pipeline into PlasmidNet SQLite.

Usage:
    # After downloading results from EC2:
    python import_mobsuite_results.py mob_results.csv amr_results.tsv

This updates the ncbi_extra table with Inc groups and mobility,
and adds AMR annotations to the amr table.
"""

import csv
import os
import sqlite3
import sys

DB_PATH = os.path.join(os.path.dirname(__file__), "data", "plsdb.db")


def import_mob_results(csv_path):
    """Import MOBsuite mob_typer results."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    # Add columns to ncbi_extra if they don't exist
    for col, coltype in [("rep_type", "TEXT"), ("predicted_mobility", "TEXT"),
                          ("relaxase_type", "TEXT"), ("mpf_type", "TEXT"),
                          ("orit_type", "TEXT"), ("mash_neighbor", "TEXT")]:
        try:
            cur.execute(f"ALTER TABLE ncbi_extra ADD COLUMN {col} {coltype} DEFAULT ''")
        except sqlite3.OperationalError:
            pass  # column exists

    updated = 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        # MOBsuite output is tab-separated
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            # Try comma-separated
            f.seek(0)
            reader = csv.DictReader(f)

        for row in reader:
            # Extract accession from sample_id or file_id
            sample = row.get("sample_id", row.get("file_id", ""))
            # The sample_id might be a file path; extract accession
            acc = os.path.basename(sample).replace(".fasta", "")

            # Map MOBsuite output columns
            rep_type = row.get("rep_type(s)", row.get("rep_type", ""))
            mobility = row.get("predicted_mobility", "")
            relaxase = row.get("relaxase_type(s)", row.get("relaxase_type", ""))
            mpf = row.get("mpf_type", "")
            orit = row.get("orit_type(s)", row.get("orit_type", ""))
            mash = row.get("mash_nearest_neighbor", "")

            # Try to match accession in ncbi_extra
            # MOBsuite uses the FASTA header which includes the full accession
            cur.execute("""
                UPDATE ncbi_extra SET
                    rep_type = ?, predicted_mobility = ?,
                    relaxase_type = ?, mpf_type = ?,
                    orit_type = ?, mash_neighbor = ?
                WHERE accession = ? OR accession LIKE ?
            """, (rep_type, mobility, relaxase, mpf, orit, mash,
                  acc, f"{acc}%"))

            if cur.rowcount > 0:
                updated += 1

    conn.commit()

    # Create indexes
    try:
        cur.execute("CREATE INDEX idx_ncbi_mobility ON ncbi_extra(predicted_mobility)")
        cur.execute("CREATE INDEX idx_ncbi_rep ON ncbi_extra(rep_type)")
    except sqlite3.OperationalError:
        pass

    conn.close()
    print(f"MOBsuite: updated {updated:,} plasmids in ncbi_extra")

    # Stats
    conn = sqlite3.connect(DB_PATH)
    typed = conn.execute(
        "SELECT COUNT(*) FROM ncbi_extra WHERE predicted_mobility != ''"
    ).fetchone()[0]
    total = conn.execute("SELECT COUNT(*) FROM ncbi_extra").fetchone()[0]
    print(f"  Typed: {typed:,} / {total:,} ({typed/total*100:.0f}%)")

    # Mobility distribution
    rows = conn.execute("""
        SELECT predicted_mobility, COUNT(*) as cnt
        FROM ncbi_extra WHERE predicted_mobility != ''
        GROUP BY predicted_mobility ORDER BY cnt DESC
    """).fetchall()
    for mob, cnt in rows:
        print(f"  {mob}: {cnt:,}")
    conn.close()


def import_amr_results(tsv_path):
    """Import AMRFinderPlus results into the amr table."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    inserted = 0
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        batch = []
        for row in reader:
            # AMRFinderPlus output columns
            acc = row.get("Contig id", row.get("Name", ""))
            # Extract accession from contig name
            if " " in acc:
                acc = acc.split()[0]

            gene = row.get("Gene symbol", "")
            gene_name = row.get("Sequence name", "")
            drug_class = row.get("Class", row.get("Subclass", ""))
            agent = row.get("Subclass", "")
            start = row.get("Start", "")
            stop = row.get("Stop", "")
            strand = row.get("Strand", "+")
            identity = row.get("% Identity to reference sequence", "")
            coverage = row.get("% Coverage of reference sequence", "")

            try:
                start_int = int(start) if start else None
                stop_int = int(stop) if stop else None
                ident_f = float(identity) if identity else None
                cov_f = float(coverage) if coverage else None
            except (ValueError, TypeError):
                start_int = stop_int = ident_f = cov_f = None

            batch.append((
                acc, gene, gene_name, drug_class, agent,
                start_int, stop_int, strand, ident_f, cov_f,
            ))
            inserted += 1

            if len(batch) >= 5000:
                cur.executemany(
                    "INSERT INTO amr VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
                batch = []

        if batch:
            cur.executemany(
                "INSERT INTO amr VALUES (?,?,?,?,?,?,?,?,?,?)", batch)

    conn.commit()
    conn.close()
    print(f"AMRFinderPlus: inserted {inserted:,} annotations into amr table")


def main():
    if len(sys.argv) < 2:
        print("Usage: python import_mobsuite_results.py <mob_results.csv> [amr_results.tsv]")
        print("\nProvide at least the MOBsuite results file.")
        sys.exit(1)

    mob_file = sys.argv[1]
    if os.path.exists(mob_file):
        print(f"Importing MOBsuite results from {mob_file}...")
        import_mob_results(mob_file)
    else:
        print(f"File not found: {mob_file}")

    if len(sys.argv) >= 3:
        amr_file = sys.argv[2]
        if os.path.exists(amr_file):
            print(f"\nImporting AMRFinderPlus results from {amr_file}...")
            import_amr_results(amr_file)
        else:
            print(f"File not found: {amr_file}")

    print("\nDone! Restart the dashboard to see updated data.")


if __name__ == "__main__":
    main()
