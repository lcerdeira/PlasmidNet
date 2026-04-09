"""
SQLite-backed data access layer for PlasmidNet dashboard.
Replaces all PLSDB API calls with local database queries.
All 72,556 plasmids — no sampling, no API rate limits.
"""

import os
import sqlite3
from pathlib import Path

DB_PATH = os.path.join(os.path.dirname(__file__), "data", "plsdb.db")

_conn = None


def get_conn():
    global _conn
    if _conn is None:
        if not os.path.exists(DB_PATH):
            raise FileNotFoundError(
                f"Database not found at {DB_PATH}. "
                "Run 'python build_db.py' to build it from PLSDB CSV files."
            )
        _conn = sqlite3.connect(DB_PATH, check_same_thread=False)
        _conn.row_factory = sqlite3.Row
        _conn.execute("PRAGMA cache_size = -64000")  # 64 MB cache
    return _conn


def q(sql, params=()):
    """Execute a query and return all rows as list of dicts."""
    cur = get_conn().execute(sql, params)
    return [dict(row) for row in cur.fetchall()]


def q1(sql, params=()):
    """Execute a query and return the first row as a dict."""
    cur = get_conn().execute(sql, params)
    row = cur.fetchone()
    return dict(row) if row else None


def scalar(sql, params=()):
    """Execute a query and return a single scalar value."""
    cur = get_conn().execute(sql, params)
    row = cur.fetchone()
    return row[0] if row else None


# ── Overview stats ──────────────────────────────────────────────────

def overview_stats():
    total = scalar("SELECT COUNT(*) FROM nuccore")
    circular = scalar("SELECT COUNT(*) FROM nuccore WHERE NUCCORE_Topology='circular'")
    linear = scalar("SELECT COUNT(*) FROM nuccore WHERE NUCCORE_Topology='linear'")
    refseq = scalar("SELECT COUNT(*) FROM nuccore WHERE LOWER(NUCCORE_Source)='refseq'")
    insdc = scalar("SELECT COUNT(*) FROM nuccore WHERE LOWER(NUCCORE_Source) IN ('insd','insdc')")
    bacteria = scalar("""
        SELECT COUNT(*) FROM nuccore n
        JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
        WHERE t.TAXONOMY_superkingdom LIKE 'Bacteria%'
    """)
    archaea = scalar("""
        SELECT COUNT(*) FROM nuccore n
        JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
        WHERE t.TAXONOMY_superkingdom LIKE 'Archaea%'
    """)
    total_amr = scalar("SELECT COUNT(*) FROM amr WHERE drug_class IS NOT NULL AND drug_class != ''")
    total_vir = scalar("SELECT COUNT(*) FROM amr WHERE drug_class IS NULL OR drug_class=''") or 0
    return {
        "total": total, "topology_circular": circular, "topology_linear": linear,
        "source_RefSeq": refseq or 0, "source_INSDC": insdc or 0,
        "kingdom_Bacteria": bacteria, "kingdom_Archaea": archaea,
        "total_annotations": total_amr + total_vir,
        "total_virulence_factors": total_vir,
        "total_bgcs": 0,
    }


# ── Taxonomy ────────────────────────────────────────────────────────

def top_genera(limit=20):
    rows = q("""
        SELECT t.TAXONOMY_genus AS genus, COUNT(*) AS cnt
        FROM nuccore n
        JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID
        WHERE t.TAXONOMY_genus != ''
        GROUP BY t.TAXONOMY_genus
        ORDER BY cnt DESC
        LIMIT ?
    """, (limit,))
    return {r["genus"]: r["cnt"] for r in rows}


# ── AMR data ────────────────────────────────────────────────────────

def amr_drug_class_counts():
    rows = q("""
        SELECT drug_class, COUNT(DISTINCT NUCCORE_ACC) AS cnt
        FROM amr
        WHERE drug_class IS NOT NULL AND drug_class != ''
        GROUP BY drug_class
        ORDER BY cnt DESC
    """)
    return {r["drug_class"]: r["cnt"] for r in rows}


def amr_gene_counts(limit=30):
    rows = q("""
        SELECT gene_symbol, COUNT(DISTINCT NUCCORE_ACC) AS cnt
        FROM amr
        WHERE gene_symbol IS NOT NULL AND gene_symbol != ''
        GROUP BY gene_symbol
        ORDER BY cnt DESC
        LIMIT ?
    """, (limit,))
    return {r["gene_symbol"]: r["cnt"] for r in rows}


# ── MOB typing / Inc groups / Mobility ──────────────────────────────

def mobility_distribution():
    rows = q("""
        SELECT predicted_mobility, COUNT(*) AS cnt
        FROM typing
        WHERE predicted_mobility != ''
        GROUP BY predicted_mobility
    """)
    return {r["predicted_mobility"]: r["cnt"] for r in rows}


def inc_group_counts(limit=20):
    """Count plasmids per Inc group (rep_type). Handles comma-separated values."""
    rows = q("SELECT rep_type FROM typing WHERE rep_type != ''")
    counts = {}
    for r in rows:
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt:
                counts[rt] = counts.get(rt, 0) + 1
    return dict(sorted(counts.items(), key=lambda x: -x[1])[:limit])


def relaxase_counts():
    rows = q("SELECT relaxase_type FROM typing WHERE relaxase_type != ''")
    counts = {}
    for r in rows:
        for rx in r["relaxase_type"].split(","):
            rx = rx.strip()
            if rx:
                counts[rx] = counts.get(rx, 0) + 1
    return dict(sorted(counts.items(), key=lambda x: -x[1]))


def mpf_counts():
    rows = q("SELECT mpf_type FROM typing WHERE mpf_type != ''")
    counts = {}
    for r in rows:
        for mp in r["mpf_type"].split(","):
            mp = mp.strip()
            if mp:
                counts[mp] = counts.get(mp, 0) + 1
    return dict(sorted(counts.items(), key=lambda x: -x[1]))


# ── Temporal / GC / Length distributions ────────────────────────────

def temporal_distribution():
    rows = q("""
        SELECT SUBSTR(NUCCORE_CreateDate, 1, 4) AS year, COUNT(*) AS cnt
        FROM nuccore
        WHERE NUCCORE_CreateDate != '' AND LENGTH(NUCCORE_CreateDate) >= 4
        GROUP BY year
        ORDER BY year
    """)
    return {r["year"]: r["cnt"] for r in rows if r["year"] and r["year"].isdigit()}


def gc_distribution():
    rows = q("""
        SELECT
            CASE
                WHEN NUCCORE_GC < 0.30 THEN '25-30'
                WHEN NUCCORE_GC < 0.35 THEN '30-35'
                WHEN NUCCORE_GC < 0.40 THEN '35-40'
                WHEN NUCCORE_GC < 0.45 THEN '40-45'
                WHEN NUCCORE_GC < 0.50 THEN '45-50'
                WHEN NUCCORE_GC < 0.55 THEN '50-55'
                WHEN NUCCORE_GC < 0.60 THEN '55-60'
                WHEN NUCCORE_GC < 0.65 THEN '60-65'
                WHEN NUCCORE_GC < 0.70 THEN '65-70'
                ELSE '70-75'
            END AS gc_range,
            COUNT(*) AS cnt
        FROM nuccore
        WHERE NUCCORE_GC > 0
        GROUP BY gc_range
        ORDER BY gc_range
    """)
    return {r["gc_range"]: r["cnt"] for r in rows}


def length_distribution():
    rows = q("""
        SELECT
            CASE
                WHEN NUCCORE_Length < 5000 THEN '0-5kb'
                WHEN NUCCORE_Length < 10000 THEN '5-10kb'
                WHEN NUCCORE_Length < 25000 THEN '10-25kb'
                WHEN NUCCORE_Length < 50000 THEN '25-50kb'
                WHEN NUCCORE_Length < 100000 THEN '50-100kb'
                WHEN NUCCORE_Length < 200000 THEN '100-200kb'
                WHEN NUCCORE_Length < 500000 THEN '200-500kb'
                ELSE '>500kb'
            END AS size_range,
            COUNT(*) AS cnt
        FROM nuccore
        WHERE NUCCORE_Length > 0
        GROUP BY size_range
        ORDER BY
            CASE size_range
                WHEN '0-5kb' THEN 1 WHEN '5-10kb' THEN 2
                WHEN '10-25kb' THEN 3 WHEN '25-50kb' THEN 4
                WHEN '50-100kb' THEN 5 WHEN '100-200kb' THEN 6
                WHEN '200-500kb' THEN 7 ELSE 8
            END
    """)
    return {r["size_range"]: r["cnt"] for r in rows}


# ── Correlations (full database, no sampling) ───────────────────────

HEAVY_METAL_CLASSES = {"ARSENIC", "COPPER", "COPPER/SILVER", "MERCURY",
                       "SILVER", "TELLURIUM", "ZINC", "LEAD", "CADMIUM",
                       "CHROMIUM"}


def build_correlations():
    """Build all correlation matrices from the full database."""
    conn = get_conn()

    # ── Pre-compute per-plasmid Inc groups + mobility ───────────────
    typing_rows = q("SELECT NUCCORE_ACC, rep_type, predicted_mobility, "
                    "PMLST_scheme, PMLST_alleles FROM typing")
    plasmid_inc = {}    # acc -> set of inc groups
    plasmid_mob = {}    # acc -> mobility string
    plasmid_pmlst = {}  # acc -> (scheme, alleles)
    for r in typing_rows:
        acc = r["NUCCORE_ACC"]
        plasmid_mob[acc] = r["predicted_mobility"] or "unknown"
        inc_set = set()
        for rt in (r["rep_type"] or "").split(","):
            rt = rt.strip()
            if rt:
                inc_set.add(rt)
        plasmid_inc[acc] = inc_set
        if r["PMLST_scheme"]:
            plasmid_pmlst[acc] = (r["PMLST_scheme"], r["PMLST_alleles"] or "")

    total_plasmids = len(plasmid_mob)

    # ── 1. AMR drug class vs Inc / Mobility ─────────────────────────
    amr_rows = q("""
        SELECT NUCCORE_ACC, drug_class, gene_symbol
        FROM amr
        WHERE drug_class IS NOT NULL AND drug_class != ''
    """)

    # Group drug classes per plasmid
    plasmid_dc = {}     # acc -> set of drug_classes
    plasmid_metals = {} # acc -> set of metal types
    for r in amr_rows:
        acc = r["NUCCORE_ACC"]
        dc = r["drug_class"].strip().upper()
        plasmid_dc.setdefault(acc, set()).add(dc)
        if dc in HEAVY_METAL_CLASSES:
            plasmid_metals.setdefault(acc, set()).add(dc)

    # AMR vs Inc
    amr_vs_inc = {}
    amr_vs_mobility = {}
    for acc, dcs in plasmid_dc.items():
        mob = plasmid_mob.get(acc, "unknown")
        incs = plasmid_inc.get(acc, set())
        for dc in dcs:
            amr_vs_inc.setdefault(dc, {})
            for inc in incs:
                amr_vs_inc[dc][inc] = amr_vs_inc[dc].get(inc, 0) + 1
            amr_vs_mobility.setdefault(dc, {})
            amr_vs_mobility[dc][mob] = amr_vs_mobility[dc].get(mob, 0) + 1

    # ── 2. Heavy metal vs Inc / Mobility ────────────────────────────
    hm_gene_counts = {}
    for r in amr_rows:
        dc = r["drug_class"].strip().upper()
        if dc in HEAVY_METAL_CLASSES:
            g = r["gene_symbol"]
            hm_gene_counts[g] = hm_gene_counts.get(g, 0) + 1

    hm_vs_inc = {}
    hm_vs_mobility = {}
    for acc, metals in plasmid_metals.items():
        mob = plasmid_mob.get(acc, "unknown")
        incs = plasmid_inc.get(acc, set())
        for metal in metals:
            hm_vs_inc.setdefault(metal, {})
            for inc in incs:
                hm_vs_inc[metal][inc] = hm_vs_inc[metal].get(inc, 0) + 1
            hm_vs_mobility.setdefault(metal, {})
            hm_vs_mobility[metal][mob] = hm_vs_mobility[metal].get(mob, 0) + 1

    # ── 3. Virulence: plasmids with virulence-associated AMR entries
    # PLSDB AMR table includes virulence factors with NULL drug_class
    vir_rows = q("""
        SELECT DISTINCT a.NUCCORE_ACC, a.gene_symbol
        FROM amr a
        WHERE (a.drug_class IS NULL OR a.drug_class = '')
          AND a.gene_symbol != ''
    """)
    vir_plasmids = {}  # acc -> list of gene names
    for r in vir_rows:
        vir_plasmids.setdefault(r["NUCCORE_ACC"], []).append(r["gene_symbol"])

    vir_gene_counts = {}
    for genes in vir_plasmids.values():
        for g in genes:
            vir_gene_counts[g] = vir_gene_counts.get(g, 0) + 1

    vir_vs_inc = {}
    vir_vs_mobility = {"conjugative": 0, "mobilizable": 0, "non-mobilizable": 0}
    for acc in vir_plasmids:
        mob = plasmid_mob.get(acc, "unknown")
        if mob in vir_vs_mobility:
            vir_vs_mobility[mob] += 1
        for inc in plasmid_inc.get(acc, set()):
            vir_vs_inc[inc] = vir_vs_inc.get(inc, 0) + 1

    # ── 4. Toxin-Antitoxin systems ─────────────────────────────────
    # Search by gene_symbol patterns + gene_name keywords across ALL
    # AMR entries (including those with empty drug_class = virulence)
    ta_sql = """
        SELECT DISTINCT NUCCORE_ACC, gene_symbol FROM amr
        WHERE LOWER(gene_symbol) IN (
            'rele','relb','higa','higb','mazf','maze','ccda','ccdb',
            'pard','pare','vapb','vapc','phd','doc','yoeb','yefm',
            'mqsr','mqsa','hipa','hipb','pema','pemk','pemi',
            'higb-1','higa-1','relb2','rele2','vapbc','ccdab'
        )
        OR LOWER(gene_symbol) LIKE 'vap%'
        OR LOWER(gene_symbol) LIKE 'ccd%'
        OR LOWER(gene_symbol) LIKE 'pem%'
        OR LOWER(gene_symbol) LIKE 'rel%'
        OR LOWER(gene_symbol) LIKE 'hig%'
        OR LOWER(gene_symbol) LIKE 'maz%'
        OR LOWER(gene_symbol) LIKE 'par%'
        OR LOWER(gene_name) LIKE '%toxin-antitoxin%'
        OR LOWER(gene_name) LIKE '%antitoxin%'
        OR LOWER(gene_name) LIKE '%addiction module%'
        OR LOWER(gene_name) LIKE '%plasmid stabilization%'
        OR LOWER(gene_name) LIKE '%post-segregational%'
    """
    ta_rows_db = q(ta_sql)
    ta_from_db = {}
    for r in ta_rows_db:
        ta_from_db.setdefault(r["NUCCORE_ACC"], []).append(r["gene_symbol"])

    ta_gene_counts = {}
    for genes in ta_from_db.values():
        for g in genes:
            ta_gene_counts[g] = ta_gene_counts.get(g, 0) + 1

    ta_vs_inc = {}
    ta_vs_mobility = {"conjugative": 0, "mobilizable": 0, "non-mobilizable": 0}
    for acc in ta_from_db:
        mob = plasmid_mob.get(acc, "unknown")
        if mob in ta_vs_mobility:
            ta_vs_mobility[mob] += 1
        for inc in plasmid_inc.get(acc, set()):
            ta_vs_inc[inc] = ta_vs_inc.get(inc, 0) + 1

    # ── 5. QAC: quaternary ammonium compound resistance ─────────────
    qac_keywords = ("QUATERNARY AMMONIUM",)
    qac_gene_prefixes = ("qac",)

    qac_plasmids = {}
    for r in amr_rows:
        dc = (r["drug_class"] or "").upper()
        g = (r["gene_symbol"] or "").lower()
        if any(kw in dc for kw in qac_keywords) or any(g.startswith(p) for p in qac_gene_prefixes):
            qac_plasmids.setdefault(r["NUCCORE_ACC"], []).append(r["gene_symbol"])

    qac_gene_counts = {}
    for genes in qac_plasmids.values():
        for g in genes:
            qac_gene_counts[g] = qac_gene_counts.get(g, 0) + 1

    qac_vs_inc = {}
    qac_vs_mobility = {"conjugative": 0, "mobilizable": 0, "non-mobilizable": 0}
    for acc in qac_plasmids:
        mob = plasmid_mob.get(acc, "unknown")
        if mob in qac_vs_mobility:
            qac_vs_mobility[mob] += 1
        for inc in plasmid_inc.get(acc, set()):
            qac_vs_inc[inc] = qac_vs_inc.get(inc, 0) + 1

    # ── 6. pMLST ───────────────────────────────────────────────────
    pmlst_schemes = {}
    pmlst_allele_freq = {}
    pmlst_vs_mobility = {}
    for acc, (scheme, alleles) in plasmid_pmlst.items():
        pmlst_schemes[scheme] = pmlst_schemes.get(scheme, 0) + 1
        mob = plasmid_mob.get(acc, "unknown")
        pmlst_vs_mobility.setdefault(scheme, {})
        pmlst_vs_mobility[scheme][mob] = pmlst_vs_mobility[scheme].get(mob, 0) + 1
        for part in alleles.split(","):
            part = part.strip()
            if part and "(-)" not in part:
                pmlst_allele_freq[part] = pmlst_allele_freq.get(part, 0) + 1

    # ── 7. Transposase elements from typing_markers ─────────────────
    # Note: phage data requires PGAP annotations (proteins.csv, 637MB)
    # which is not included. We report MOB-associated element counts.
    mob_element_counts = {}
    mob_per_plasmid = {}
    tm_rows = q("SELECT NUCCORE_ACC, element FROM typing_markers")
    for r in tm_rows:
        el = r["element"]
        mob_element_counts[el] = mob_element_counts.get(el, 0) + 1
        mob_per_plasmid.setdefault(r["NUCCORE_ACC"], set()).add(el)

    phage_bins = {}  # empty = not available
    phage_mob = {"conjugative": [], "mobilizable": [], "non-mobilizable": []}
    for acc, elements in mob_per_plasmid.items():
        mob = plasmid_mob.get(acc, "unknown")
        if mob in phage_mob:
            phage_mob[mob].append(len(elements))

    phage_mob_summary = {}
    for mob, vals in phage_mob.items():
        if vals:
            phage_mob_summary[mob] = {
                "mean_phage_elements": round(sum(vals) / len(vals), 2),
                "count": len(vals), "total_elements": sum(vals),
            }

    return {
        "sample_size": total_plasmids,
        "amr_vs_inc": amr_vs_inc,
        "amr_vs_mobility": amr_vs_mobility,
        "heavy_metal_genes": hm_gene_counts,
        "heavy_metal_vs_inc": hm_vs_inc,
        "heavy_metal_vs_mobility": hm_vs_mobility,
        "phage_distribution": phage_bins,
        "transposase_distribution": {},
        "phage_vs_mobility": phage_mob_summary,
        "pmlst_schemes": pmlst_schemes,
        "pmlst_alleles": pmlst_allele_freq,
        "pmlst_vs_mobility": pmlst_vs_mobility,
        "vir_total": len(vir_plasmids),
        "vir_gene_counts": vir_gene_counts,
        "vir_vs_inc": vir_vs_inc,
        "vir_vs_mobility": vir_vs_mobility,
        "ta_total": len(ta_from_db),
        "ta_gene_counts": ta_gene_counts,
        "ta_vs_inc": ta_vs_inc,
        "ta_vs_mobility": ta_vs_mobility,
        "qac_total": len(qac_plasmids),
        "qac_gene_counts": qac_gene_counts,
        "qac_vs_inc": qac_vs_inc,
        "qac_vs_mobility": qac_vs_mobility,
    }


# ── Plasmid viewer lookup ──────────────────────────────────────────

def plasmid_summary(accession):
    """Get full plasmid data for the viewer (replaces PLSDB API /summary)."""
    nuc = q1("SELECT * FROM nuccore WHERE NUCCORE_ACC = ?", (accession,))
    if not nuc:
        return None

    tax = q1("SELECT * FROM taxonomy WHERE TAXONOMY_UID = ?",
             (nuc.get("TAXONOMY_UID"),)) if nuc.get("TAXONOMY_UID") else {}

    typ = q1("SELECT * FROM typing WHERE NUCCORE_ACC = ?", (accession,))

    amr_feats = q("""
        SELECT gene_symbol, gene_name, drug_class, antimicrobial_agent,
               input_gene_start, input_gene_stop, strand_orientation,
               sequence_identity, coverage_percentage
        FROM amr WHERE NUCCORE_ACC = ?
    """, (accession,))

    mob_feats = q("""
        SELECT element, sstart, send, sstrand, pident, qcovhsp, biomarker
        FROM typing_markers WHERE NUCCORE_ACC = ?
    """, (accession,))

    return {
        "nuccore": nuc,
        "taxonomy": tax or {},
        "typing": typ or {},
        "amr": amr_feats,
        "mob": mob_feats,
    }
