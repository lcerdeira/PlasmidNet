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

def has_ncbi_extra():
    """Check if ncbi_extra table exists."""
    try:
        scalar("SELECT COUNT(*) FROM ncbi_extra")
        return True
    except Exception:
        return False


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
        "total_amr": total_amr,
        "total_bgcs": 0,
        "ncbi_extra": scalar("SELECT COUNT(*) FROM ncbi_extra") if has_ncbi_extra() else 0,
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


def amr_genes_with_class(limit=30):
    """Return top AMR genes with their drug class from the database."""
    rows = q("""
        SELECT gene_symbol, drug_class, COUNT(DISTINCT NUCCORE_ACC) AS cnt
        FROM amr
        WHERE gene_symbol IS NOT NULL AND gene_symbol != ''
          AND drug_class IS NOT NULL AND drug_class != ''
        GROUP BY gene_symbol, drug_class
        ORDER BY cnt DESC
        LIMIT ?
    """, (limit * 2,))
    # Deduplicate genes (keep highest count class)
    seen = {}
    for r in rows:
        g = r["gene_symbol"]
        if g not in seen:
            seen[g] = {"gene": g, "drug_class": r["drug_class"], "count": r["cnt"]}
    result = sorted(seen.values(), key=lambda x: -x["count"])[:limit]
    return result


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

    # ── 4. Toxin-Antitoxin systems (from pgap_features table) ──────
    ta_rows_db = q("""
        SELECT NUCCORE_ACC, gene FROM pgap_features
        WHERE category = 'TA'
    """)
    ta_from_db = {}
    for r in ta_rows_db:
        ta_from_db.setdefault(r["NUCCORE_ACC"], []).append(r["gene"])

    ta_gene_counts = {}
    for genes in ta_from_db.values():
        for g in genes:
            if g:
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

    # ── 7. Phage & transposase (from pgap_features table) ───────────
    phage_per_plasmid = {}
    trans_per_plasmid = {}
    pgap_rows = q("""
        SELECT NUCCORE_ACC, category, COUNT(*) as cnt
        FROM pgap_features
        WHERE category IN ('PHAGE', 'TRANSPOSASE')
        GROUP BY NUCCORE_ACC, category
    """)
    for r in pgap_rows:
        if r["category"] == "PHAGE":
            phage_per_plasmid[r["NUCCORE_ACC"]] = r["cnt"]
        else:
            trans_per_plasmid[r["NUCCORE_ACC"]] = r["cnt"]

    phage_bins = {"0": 0, "1-2": 0, "3-5": 0, "6-10": 0, ">10": 0}
    phage_bins["0"] = total_plasmids - len(phage_per_plasmid)
    for cnt in phage_per_plasmid.values():
        if cnt <= 2: phage_bins["1-2"] += 1
        elif cnt <= 5: phage_bins["3-5"] += 1
        elif cnt <= 10: phage_bins["6-10"] += 1
        else: phage_bins[">10"] += 1

    trans_bins = {"0": 0, "1-3": 0, "4-8": 0, "9-15": 0, ">15": 0}
    trans_bins["0"] = total_plasmids - len(trans_per_plasmid)
    for cnt in trans_per_plasmid.values():
        if cnt <= 3: trans_bins["1-3"] += 1
        elif cnt <= 8: trans_bins["4-8"] += 1
        elif cnt <= 15: trans_bins["9-15"] += 1
        else: trans_bins[">15"] += 1

    phage_mob = {"conjugative": [], "mobilizable": [], "non-mobilizable": []}
    for acc, cnt in phage_per_plasmid.items():
        mob = plasmid_mob.get(acc, "unknown")
        if mob in phage_mob:
            phage_mob[mob].append(cnt)

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
        "transposase_distribution": trans_bins,
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


# ── Geography queries ──────────────────────────────────────────────

def geo_plasmid_data():
    """
    Return per-plasmid geographic data joined with typing and AMR info.
    Used for global maps, country comparisons, and temporal animations.
    """
    rows = q("""
        SELECT
            n.NUCCORE_ACC, n.NUCCORE_CreateDate,
            SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year,
            b.lat, b.lng, b.country,
            t.rep_type, t.predicted_mobility, t.mpf_type, t.relaxase_type
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country != '' AND b.lat IS NOT NULL
    """)
    return rows


def geo_country_summary():
    """Plasmid counts per country with Inc/mobility breakdown."""
    rows = q("""
        SELECT b.country, t.predicted_mobility, COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country != ''
        GROUP BY b.country, t.predicted_mobility
        ORDER BY cnt DESC
    """)
    return rows


def geo_country_inc_top(limit=15):
    """Top Inc groups per country (top countries only)."""
    rows = q("""
        SELECT b.country, t.rep_type, COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country != '' AND t.rep_type != ''
        GROUP BY b.country
        ORDER BY cnt DESC
        LIMIT ?
    """, (limit,))
    # Get top countries
    top_countries = [r["country"] for r in rows]

    # Now get Inc distribution for those countries
    if not top_countries:
        return {}
    placeholders = ",".join("?" * len(top_countries))
    detail = q(f"""
        SELECT b.country, t.rep_type
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country IN ({placeholders}) AND t.rep_type != ''
    """, tuple(top_countries))

    # Parse comma-separated rep_types
    result = {}
    for r in detail:
        country = r["country"]
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt:
                result.setdefault(country, {})
                result[country][rt] = result[country].get(rt, 0) + 1
    return result


def geo_temporal_data():
    """Year-by-year plasmid appearance per country for animation."""
    rows = q("""
        SELECT
            SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year,
            b.country, b.lat, b.lng,
            COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        WHERE b.country != '' AND b.lat IS NOT NULL
          AND LENGTH(n.NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2000
        GROUP BY year, b.country
        ORDER BY year, cnt DESC
    """)
    return rows


def geo_temporal_inc(inc_group="IncFIB"):
    """Year-by-year spread of a specific Inc group across countries."""
    rows = q("""
        SELECT
            SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year,
            b.country, b.lat, b.lng,
            COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country != '' AND b.lat IS NOT NULL
          AND t.rep_type LIKE ?
          AND LENGTH(n.NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2000
        GROUP BY year, b.country
        ORDER BY year, cnt DESC
    """, (f"%{inc_group}%",))
    return rows


def geo_temporal_mobility(mobility_type="conjugative"):
    """Year-by-year spread of a mobility type across countries."""
    rows = q("""
        SELECT
            SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year,
            b.country, b.lat, b.lng,
            COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country != '' AND b.lat IS NOT NULL
          AND t.predicted_mobility = ?
          AND LENGTH(n.NUCCORE_CreateDate) >= 4
          AND CAST(SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS INTEGER) >= 2000
        GROUP BY year, b.country
        ORDER BY year, cnt DESC
    """, (mobility_type,))
    return rows


# ── Host / source correlations ─────────────────────────────────────

def host_category_counts():
    """Count plasmids per host category."""
    rows = q("""
        SELECT b.host_category, COUNT(DISTINCT n.NUCCORE_ACC) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        WHERE b.host_category != '' AND b.host_category != 'Other'
        GROUP BY b.host_category
        ORDER BY cnt DESC
    """)
    return {r["host_category"]: r["cnt"] for r in rows}


def host_vs_mobility():
    """Mobility distribution per host category."""
    rows = q("""
        SELECT b.host_category, t.predicted_mobility, COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.host_category != '' AND b.host_category != 'Other'
        GROUP BY b.host_category, t.predicted_mobility
    """)
    result = {}
    for r in rows:
        result.setdefault(r["host_category"], {})
        result[r["host_category"]][r["predicted_mobility"]] = r["cnt"]
    return result


def host_vs_inc(limit=10):
    """Top Inc groups per host category."""
    rows = q("""
        SELECT b.host_category, t.rep_type
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.host_category != '' AND b.host_category != 'Other'
          AND t.rep_type != ''
    """)
    result = {}
    for r in rows:
        cat = r["host_category"]
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt:
                result.setdefault(cat, {})
                result[cat][rt] = result[cat].get(rt, 0) + 1
    return result


def host_vs_amr_class():
    """Top AMR drug classes per host category."""
    rows = q("""
        SELECT b.host_category, a.drug_class, COUNT(DISTINCT n.NUCCORE_ACC) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN amr a ON n.NUCCORE_ACC = a.NUCCORE_ACC
        WHERE b.host_category != '' AND b.host_category != 'Other'
          AND a.drug_class IS NOT NULL AND a.drug_class != ''
        GROUP BY b.host_category, a.drug_class
        ORDER BY cnt DESC
    """)
    result = {}
    for r in rows:
        result.setdefault(r["host_category"], {})
        result[r["host_category"]][r["drug_class"]] = r["cnt"]
    return result


# ── Analytics queries ──────────────────────────────────────────────

TOP_COUNTRIES = ("China", "USA", "United Kingdom", "South Korea", "Japan",
                 "Germany", "Australia", "Canada", "France", "Spain")
TOP_GENERA = ("Escherichia", "Klebsiella", "Salmonella", "Staphylococcus",
              "Acinetobacter", "Enterococcus")


def analytics_matched_comparison():
    """
    Matched comparison: mobility ratios within same species across countries.
    Controls for species confound.
    """
    ph = ",".join("?" * len(TOP_COUNTRIES))
    pg = ",".join("?" * len(TOP_GENERA))
    rows = q(f"""
        SELECT tx.TAXONOMY_genus AS genus, b.country,
               t.predicted_mobility AS mobility, COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        JOIN taxonomy tx ON n.TAXONOMY_UID = tx.TAXONOMY_UID
        WHERE b.country IN ({ph})
          AND tx.TAXONOMY_genus IN ({pg})
          AND t.predicted_mobility != ''
        GROUP BY tx.TAXONOMY_genus, b.country, t.predicted_mobility
    """, TOP_COUNTRIES + TOP_GENERA)
    return rows


def analytics_rarefaction(n_iter=50, sample_sizes=None):
    """
    Rarefaction: subsample each country to equal n and compute
    conjugative fraction. Returns data for rarefaction curves.
    """
    import random
    random.seed(42)

    if sample_sizes is None:
        sample_sizes = [50, 100, 200, 500, 1000, 2000]

    ph = ",".join("?" * len(TOP_COUNTRIES))
    rows = q(f"""
        SELECT b.country, t.predicted_mobility
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE b.country IN ({ph}) AND t.predicted_mobility != ''
    """, TOP_COUNTRIES)

    # Group by country
    country_data = {}
    for r in rows:
        country_data.setdefault(r["country"], []).append(r["predicted_mobility"])

    results = []
    for country, mobilities in country_data.items():
        total = len(mobilities)
        for n in sample_sizes:
            if n > total:
                continue
            fractions = []
            for _ in range(n_iter):
                sample = random.sample(mobilities, n)
                conj = sum(1 for m in sample if m == "conjugative")
                fractions.append(conj / n)
            mean_frac = sum(fractions) / len(fractions)
            std_frac = (sum((f - mean_frac) ** 2 for f in fractions) / len(fractions)) ** 0.5
            results.append({
                "country": country, "n": n, "total": total,
                "conj_mean": round(mean_frac, 4),
                "conj_std": round(std_frac, 4),
            })
    return results


def analytics_feature_matrix():
    """
    Build feature matrix for ML: one row per plasmid with
    country, genus, host_category, year, length, gc, mobility.
    """
    rows = q("""
        SELECT n.NUCCORE_ACC,
               b.country, b.host_category,
               tx.TAXONOMY_genus AS genus,
               SUBSTR(n.NUCCORE_CreateDate, 1, 4) AS year,
               n.NUCCORE_Length AS length,
               n.NUCCORE_GC AS gc,
               t.predicted_mobility AS mobility,
               t.rep_type
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        JOIN taxonomy tx ON n.TAXONOMY_UID = tx.TAXONOMY_UID
        WHERE b.country != '' AND t.predicted_mobility != ''
          AND tx.TAXONOMY_genus != ''
          AND LENGTH(n.NUCCORE_CreateDate) >= 4
          AND n.NUCCORE_Length > 0
    """)
    return rows


def analytics_simpson_paradox():
    """
    Detect Simpson's paradox: cases where the mobility ratio for an
    Inc group flips direction when stratifying by species.

    Returns list of {inc_group, overall_conj_pct, species_conj_pcts}.
    """
    ph = ",".join("?" * len(TOP_COUNTRIES))
    rows = q(f"""
        SELECT t.rep_type, tx.TAXONOMY_genus AS genus,
               t.predicted_mobility AS mobility, COUNT(*) AS cnt
        FROM nuccore n
        JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        JOIN taxonomy tx ON n.TAXONOMY_UID = tx.TAXONOMY_UID
        WHERE b.country IN ({ph}) AND t.rep_type != ''
          AND t.predicted_mobility IN ('conjugative', 'non-mobilizable')
          AND tx.TAXONOMY_genus IN ('Escherichia', 'Klebsiella', 'Salmonella',
                                    'Staphylococcus', 'Acinetobacter', 'Enterococcus')
        GROUP BY t.rep_type, tx.TAXONOMY_genus, t.predicted_mobility
    """, TOP_COUNTRIES)

    # Parse comma-separated rep_types and aggregate
    inc_genus_mob = {}  # (inc, genus) -> {conj: n, non: n}
    for r in rows:
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if not rt:
                continue
            key = (rt, r["genus"])
            inc_genus_mob.setdefault(key, {"conjugative": 0, "non-mobilizable": 0})
            inc_genus_mob[key][r["mobility"]] += r["cnt"]

    # Find Inc groups where overall % conj differs from species-level
    inc_overall = {}
    for (inc, genus), mobs in inc_genus_mob.items():
        inc_overall.setdefault(inc, {"conjugative": 0, "non-mobilizable": 0})
        inc_overall[inc]["conjugative"] += mobs["conjugative"]
        inc_overall[inc]["non-mobilizable"] += mobs["non-mobilizable"]

    paradoxes = []
    for inc, overall in inc_overall.items():
        total = overall["conjugative"] + overall["non-mobilizable"]
        if total < 50:
            continue
        overall_pct = overall["conjugative"] / total

        # Per-species
        species_pcts = {}
        for (i, genus), mobs in inc_genus_mob.items():
            if i != inc:
                continue
            sp_total = mobs["conjugative"] + mobs["non-mobilizable"]
            if sp_total < 10:
                continue
            species_pcts[genus] = mobs["conjugative"] / sp_total

        if not species_pcts:
            continue

        # Check for paradox or large species-level divergence
        above = sum(1 for p in species_pcts.values() if p > 0.5)
        below = len(species_pcts) - above
        is_paradox = ((overall_pct > 0.5 and below > above) or
                      (overall_pct < 0.5 and above > below))

        # Also flag if species differ by >20 percentage points
        vals = list(species_pcts.values())
        max_spread = max(vals) - min(vals) if len(vals) >= 2 else 0

        if is_paradox or max_spread > 0.20:
            paradoxes.append({
                "inc_group": inc,
                "overall_conj_pct": round(overall_pct * 100, 1),
                "total": total,
                "species": {k: round(v * 100, 1) for k, v in
                            sorted(species_pcts.items(), key=lambda x: -x[1])},
                "paradox": is_paradox,
                "spread": round(max_spread * 100, 1),
            })

    return sorted(paradoxes, key=lambda x: -x["spread"])[:15]


# ── PubMed / network / clustering queries ──────────────────────────

def pubmed_stats():
    """Count plasmids with associated PMIDs and top cited PMIDs."""
    try:
        total = scalar("SELECT COUNT(*) FROM typing WHERE associated_pmid != ''") or 0
        rows = q("SELECT associated_pmid FROM typing WHERE associated_pmid != ''")
    except Exception:
        return {"total_with_pmid": 0, "top_pmids": []}
    pmid_counts = {}
    for r in rows:
        for pmid in r["associated_pmid"].split(";"):
            pmid = pmid.strip()
            if pmid:
                pmid_counts[pmid] = pmid_counts.get(pmid, 0) + 1
    top = sorted(pmid_counts.items(), key=lambda x: -x[1])[:20]
    return {"total_with_pmid": total, "top_pmids": top}


def integron_analysis():
    """
    Detect class 1 integrons and analyse gene cassette composition.
    Class 1 signature: qacEdelta1 + sul1 (3' conserved segment).
    """
    # Plasmids with integron markers
    class1_total = scalar("""
        SELECT COUNT(DISTINCT a1.NUCCORE_ACC) FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        WHERE a1.gene_symbol = 'qacEdelta1' AND a2.gene_symbol = 'sul1'
    """) or 0

    qac_total = scalar(
        "SELECT COUNT(DISTINCT NUCCORE_ACC) FROM amr WHERE gene_symbol='qacEdelta1'"
    ) or 0

    # Gene cassettes within 5kb of qacEdelta1
    cassette_genes = q("""
        SELECT a2.gene_symbol AS gene, a2.drug_class AS drug_class,
               COUNT(DISTINCT a1.NUCCORE_ACC) AS cnt
        FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        WHERE a1.gene_symbol = 'qacEdelta1' AND a2.gene_symbol != 'qacEdelta1'
          AND a2.gene_symbol != '' AND a1.input_gene_start IS NOT NULL
          AND a2.input_gene_start IS NOT NULL
          AND ABS(CAST(a1.input_gene_start AS INT) - CAST(a2.input_gene_start AS INT)) < 5000
        GROUP BY a2.gene_symbol ORDER BY cnt DESC LIMIT 20
    """)

    # Integron vs mobility
    integron_mob = q("""
        SELECT t.predicted_mobility, COUNT(DISTINCT a1.NUCCORE_ACC) AS cnt
        FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        JOIN typing t ON a1.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE a1.gene_symbol = 'qacEdelta1' AND a2.gene_symbol = 'sul1'
          AND t.predicted_mobility != ''
        GROUP BY t.predicted_mobility
    """)
    mob_dist = {r["predicted_mobility"]: r["cnt"] for r in integron_mob}

    # Integron vs Inc group
    integron_inc = q("""
        SELECT t.rep_type, COUNT(DISTINCT a1.NUCCORE_ACC) AS cnt
        FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        JOIN typing t ON a1.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE a1.gene_symbol = 'qacEdelta1' AND a2.gene_symbol = 'sul1'
          AND t.rep_type != ''
        GROUP BY t.rep_type ORDER BY cnt DESC
    """)
    # Parse comma-separated rep_types
    inc_counts = {}
    for r in integron_inc:
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt:
                inc_counts[rt] = inc_counts.get(rt, 0) + r["cnt"]

    return {
        "class1_total": class1_total,
        "qac_total": qac_total,
        "cassette_genes": cassette_genes,
        "mobility": mob_dist,
        "inc_groups": dict(sorted(inc_counts.items(), key=lambda x: -x[1])[:15]),
    }


def comobilization_data():
    """
    Get data for co-mobilization analysis: mobilizable plasmids
    that share a host (same BIOSAMPLE_UID) with conjugative plasmids.
    """
    # Which Inc groups of conjugative plasmids co-exist with mobilizable ones
    rows = q("""
        SELECT t_conj.rep_type AS conj_inc, t_mob.rep_type AS mob_inc,
               COUNT(*) AS cnt
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t_conj ON n1.NUCCORE_ACC = t_conj.NUCCORE_ACC
        JOIN typing t_mob ON n2.NUCCORE_ACC = t_mob.NUCCORE_ACC
        WHERE t_conj.predicted_mobility = 'conjugative'
          AND t_mob.predicted_mobility = 'mobilizable'
          AND n1.BIOSAMPLE_UID > 0
          AND t_conj.rep_type != '' AND t_mob.rep_type != ''
        GROUP BY t_conj.rep_type, t_mob.rep_type
        ORDER BY cnt DESC
        LIMIT 200
    """)

    # Parse and aggregate
    pairs = {}
    for r in rows:
        for c_inc in r["conj_inc"].split(","):
            c_inc = c_inc.strip()
            if not c_inc:
                continue
            for m_inc in r["mob_inc"].split(","):
                m_inc = m_inc.strip()
                if not m_inc:
                    continue
                key = (c_inc, m_inc)
                pairs[key] = pairs.get(key, 0) + r["cnt"]

    # Top pairs
    top_pairs = sorted(pairs.items(), key=lambda x: -x[1])[:30]

    # Summary: mobilizable plasmids in same host as conjugative
    mob_with_conj = scalar("""
        SELECT COUNT(DISTINCT n2.NUCCORE_ACC)
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t1 ON n1.NUCCORE_ACC = t1.NUCCORE_ACC
        JOIN typing t2 ON n2.NUCCORE_ACC = t2.NUCCORE_ACC
        WHERE t1.predicted_mobility = 'conjugative'
          AND t2.predicted_mobility = 'mobilizable'
          AND n1.BIOSAMPLE_UID > 0
    """) or 0

    total_mob = scalar(
        "SELECT COUNT(*) FROM typing WHERE predicted_mobility='mobilizable'"
    ) or 1

    return {
        "mob_with_conj": mob_with_conj,
        "total_mobilizable": total_mob,
        "pct_colocated": round(mob_with_conj / total_mob * 100, 1),
        "top_pairs": [(c, m, n) for (c, m), n in top_pairs],
    }


def comobilization_ml_data():
    """
    Build feature matrix for co-mobilization ML:
    Predict whether a mobilizable plasmid is co-located with a conjugative one
    based on relaxase type, Inc group, host genus, length, etc.
    """
    # Positive examples: mobilizable plasmids that ARE co-located with conjugative
    positives = q("""
        SELECT DISTINCT n_mob.NUCCORE_ACC AS acc,
               t_mob.relaxase_type, t_mob.rep_type AS mob_inc,
               t_mob.orit_type,
               n_mob.NUCCORE_Length AS length, n_mob.NUCCORE_GC AS gc,
               tx.TAXONOMY_genus AS genus,
               t_conj.mpf_type AS conj_mpf, t_conj.rep_type AS conj_inc
        FROM nuccore n_mob
        JOIN nuccore n_conj ON n_mob.BIOSAMPLE_UID = n_conj.BIOSAMPLE_UID
                           AND n_mob.NUCCORE_ACC != n_conj.NUCCORE_ACC
        JOIN typing t_mob ON n_mob.NUCCORE_ACC = t_mob.NUCCORE_ACC
        JOIN typing t_conj ON n_conj.NUCCORE_ACC = t_conj.NUCCORE_ACC
        JOIN taxonomy tx ON n_mob.TAXONOMY_UID = tx.TAXONOMY_UID
        WHERE t_mob.predicted_mobility = 'mobilizable'
          AND t_conj.predicted_mobility = 'conjugative'
          AND n_mob.BIOSAMPLE_UID > 0
          AND t_mob.relaxase_type != ''
    """)

    # Negative examples: mobilizable plasmids NOT co-located with conjugative
    negatives = q("""
        SELECT t_mob.NUCCORE_ACC AS acc,
               t_mob.relaxase_type, t_mob.rep_type AS mob_inc,
               t_mob.orit_type,
               n.NUCCORE_Length AS length, n.NUCCORE_GC AS gc,
               tx.TAXONOMY_genus AS genus
        FROM typing t_mob
        JOIN nuccore n ON t_mob.NUCCORE_ACC = n.NUCCORE_ACC
        JOIN taxonomy tx ON n.TAXONOMY_UID = tx.TAXONOMY_UID
        WHERE t_mob.predicted_mobility = 'mobilizable'
          AND t_mob.relaxase_type != ''
          AND t_mob.NUCCORE_ACC NOT IN (
              SELECT DISTINCT n2.NUCCORE_ACC
              FROM nuccore n2
              JOIN nuccore n3 ON n2.BIOSAMPLE_UID = n3.BIOSAMPLE_UID
                             AND n2.NUCCORE_ACC != n3.NUCCORE_ACC
              JOIN typing t3 ON n3.NUCCORE_ACC = t3.NUCCORE_ACC
              WHERE t3.predicted_mobility = 'conjugative'
                AND n2.BIOSAMPLE_UID > 0
          )
        LIMIT 10000
    """)

    return positives, negatives


def relaxase_t4ss_compatibility():
    """
    Observed relaxase x T4SS co-occurrence vs expected (from literature).
    Returns data for a comparison chart.
    """
    # Known rules
    KNOWN = {
        ("MOBF", "MPF_F"): True,
        ("MOBP", "MPF_T"): True, ("MOBP", "MPF_F"): True,
        ("MOBQ", "MPF_T"): True, ("MOBQ", "MPF_F"): True,
        ("MOBH", "MPF_I"): True,
        ("MOBC", "MPF_T"): True,
        ("MOBV", "MPF_T"): True,
    }

    # Observed
    rows = q("""
        SELECT t_conj.mpf_type AS mpf, t_mob.relaxase_type AS relaxase,
               COUNT(*) AS cnt
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t_conj ON n1.NUCCORE_ACC = t_conj.NUCCORE_ACC
        JOIN typing t_mob ON n2.NUCCORE_ACC = t_mob.NUCCORE_ACC
        WHERE t_conj.predicted_mobility = 'conjugative'
          AND t_mob.predicted_mobility = 'mobilizable'
          AND n1.BIOSAMPLE_UID > 0
          AND t_conj.mpf_type != '' AND t_mob.relaxase_type != ''
        GROUP BY mpf, relaxase
    """)

    # Parse comma-separated relaxase types
    observed = {}
    for r in rows:
        mpf = r["mpf"].strip()
        for relax in r["relaxase"].split(","):
            relax = relax.strip()
            if relax and mpf:
                key = (relax, mpf)
                observed[key] = observed.get(key, 0) + r["cnt"]

    result = []
    all_relax = sorted(set(k[0] for k in observed))
    all_mpf = sorted(set(k[1] for k in observed))

    for relax in all_relax:
        for mpf in all_mpf:
            cnt = observed.get((relax, mpf), 0)
            known = KNOWN.get((relax, mpf), False)
            result.append({
                "relaxase": relax, "mpf": mpf, "count": cnt,
                "known_compatible": known,
            })

    return result, all_relax, all_mpf


def amr_cooccurrence(min_count=50):
    """
    AMR gene co-occurrence: which resistance genes appear together
    on the same plasmid. Returns list of (gene1, gene2, count).
    """
    rows = q("""
        SELECT a1.gene_symbol AS gene1, a2.gene_symbol AS gene2,
               COUNT(DISTINCT a1.NUCCORE_ACC) AS cnt
        FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        WHERE a1.gene_symbol < a2.gene_symbol
          AND a1.drug_class != '' AND a2.drug_class != ''
          AND a1.gene_symbol != '' AND a2.gene_symbol != ''
        GROUP BY a1.gene_symbol, a2.gene_symbol
        HAVING cnt >= ?
        ORDER BY cnt DESC
        LIMIT 100
    """, (min_count,))
    return rows


def host_species_sharing():
    """
    Which host species share the same plasmid Inc groups.
    Returns co-occurrence data for a network.
    """
    rows = q("""
        SELECT t1.TAXONOMY_genus AS genus1, t2.TAXONOMY_genus AS genus2,
               ty.rep_type, COUNT(*) AS cnt
        FROM nuccore n1
        JOIN taxonomy t1 ON n1.TAXONOMY_UID = t1.TAXONOMY_UID
        JOIN typing ty ON n1.NUCCORE_ACC = ty.NUCCORE_ACC
        JOIN nuccore n2 ON n2.TAXONOMY_UID != n1.TAXONOMY_UID
        JOIN taxonomy t2 ON n2.TAXONOMY_UID = t2.TAXONOMY_UID
        JOIN typing ty2 ON n2.NUCCORE_ACC = ty2.NUCCORE_ACC
        WHERE t1.TAXONOMY_genus < t2.TAXONOMY_genus
          AND ty.rep_type = ty2.rep_type
          AND ty.rep_type != ''
          AND t1.TAXONOMY_genus IN ('Escherichia','Klebsiella','Salmonella',
                                    'Enterococcus','Staphylococcus','Acinetobacter')
          AND t2.TAXONOMY_genus IN ('Escherichia','Klebsiella','Salmonella',
                                    'Enterococcus','Staphylococcus','Acinetobacter')
        GROUP BY t1.TAXONOMY_genus, t2.TAXONOMY_genus
        ORDER BY cnt DESC
        LIMIT 30
    """)
    return rows


# ── IS / Transposon families ───────────────────────────────────────

def is_family_counts():
    """Count IS element families from PGAP product descriptions."""
    import re
    rows = q("""
        SELECT product, COUNT(*) as cnt, COUNT(DISTINCT NUCCORE_ACC) as plasmids
        FROM pgap_features
        WHERE category = 'TRANSPOSASE'
        GROUP BY product ORDER BY cnt DESC
    """)
    # Parse IS family from product name
    families = {}
    for r in rows:
        product = r["product"]
        # Extract IS family: "IS6-like element IS26 family transposase" -> "IS26"
        # "IS3 family transposase" -> "IS3"
        m = re.search(r"(IS\w+)\s+family", product)
        if m:
            fam = m.group(1)
        elif "integrase" in product.lower() or "intI" in product.lower():
            fam = "Integrase"
        elif "recombinase" in product.lower():
            fam = "Recombinase"
        elif "Tn3" in product:
            fam = "Tn3"
        elif "Rpn" in product:
            fam = "Rpn"
        else:
            fam = "Other"
        families.setdefault(fam, {"count": 0, "plasmids": 0})
        families[fam]["count"] += r["cnt"]
        families[fam]["plasmids"] += r["plasmids"]
    return dict(sorted(families.items(), key=lambda x: -x[1]["count"]))


def is_family_by_mobility():
    """IS family distribution by plasmid mobility."""
    import re
    rows = q("""
        SELECT p.product, t.predicted_mobility, COUNT(*) as cnt
        FROM pgap_features p
        JOIN typing t ON p.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE p.category = 'TRANSPOSASE' AND t.predicted_mobility != ''
        GROUP BY p.product, t.predicted_mobility
    """)
    families_mob = {}
    for r in rows:
        product = r["product"]
        m = re.search(r"(IS\w+)\s+family", product)
        if m:
            fam = m.group(1)
        elif "Tn3" in product:
            fam = "Tn3"
        elif "integrase" in product.lower():
            fam = "Integrase"
        else:
            continue
        families_mob.setdefault(fam, {})
        families_mob[fam][r["predicted_mobility"]] = \
            families_mob[fam].get(r["predicted_mobility"], 0) + r["cnt"]
    return families_mob


def is_family_by_inc():
    """Top IS families per Inc group."""
    import re
    rows = q("""
        SELECT p.product, t.rep_type, COUNT(*) as cnt
        FROM pgap_features p
        JOIN typing t ON p.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE p.category = 'TRANSPOSASE' AND t.rep_type != ''
        GROUP BY p.product, t.rep_type
        ORDER BY cnt DESC
        LIMIT 5000
    """)
    # Parse and aggregate
    inc_is = {}
    for r in rows:
        m = re.search(r"(IS\w+)\s+family", r["product"])
        fam = m.group(1) if m else None
        if not fam:
            continue
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt and "IncF" in rt or rt in ("IncI1", "IncN", "IncHI2", "IncC",
                                              "ColRNAI", "IncR", "IncX4"):
                inc_is.setdefault(rt, {})
                inc_is[rt][fam] = inc_is[rt].get(fam, 0) + r["cnt"]
    return inc_is


# ── Transposon / Integron / KPC analytics ──────────────────────────

def kpc_context_analysis():
    """Analyse the genetic context around blaKPC genes database-wide."""
    # KPC variants
    variants = q("""
        SELECT gene_symbol, COUNT(DISTINCT NUCCORE_ACC) as cnt
        FROM amr WHERE gene_symbol LIKE '%KPC%'
        GROUP BY gene_symbol ORDER BY cnt DESC LIMIT 20
    """)

    # Gene context (within 5kb)
    context = q("""
        SELECT a2.gene_symbol AS gene, a2.drug_class,
               COUNT(DISTINCT a1.NUCCORE_ACC) AS cnt
        FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        WHERE a1.gene_symbol LIKE '%KPC%' AND a2.gene_symbol NOT LIKE '%KPC%'
          AND a2.gene_symbol != ''
          AND a1.input_gene_start IS NOT NULL AND a2.input_gene_start IS NOT NULL
          AND ABS(CAST(a1.input_gene_start AS INT) - CAST(a2.input_gene_start AS INT)) < 5000
        GROUP BY a2.gene_symbol ORDER BY cnt DESC LIMIT 20
    """)

    # KPC by Inc group
    kpc_inc = {}
    rows = q("""
        SELECT t.rep_type, COUNT(DISTINCT a.NUCCORE_ACC) as cnt
        FROM amr a JOIN typing t ON a.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE a.gene_symbol LIKE '%KPC%' AND t.rep_type != ''
        GROUP BY t.rep_type ORDER BY cnt DESC
    """)
    for r in rows:
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt:
                kpc_inc[rt] = kpc_inc.get(rt, 0) + r["cnt"]

    # KPC by mobility
    kpc_mob = q("""
        SELECT t.predicted_mobility, COUNT(DISTINCT a.NUCCORE_ACC) as cnt
        FROM amr a JOIN typing t ON a.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE a.gene_symbol LIKE '%KPC%'
        GROUP BY t.predicted_mobility ORDER BY cnt DESC
    """)

    total_kpc = scalar("SELECT COUNT(DISTINCT NUCCORE_ACC) FROM amr WHERE gene_symbol LIKE '%KPC%'")

    return {
        "total": total_kpc,
        "variants": variants,
        "context": context,
        "inc_groups": dict(sorted(kpc_inc.items(), key=lambda x: -x[1])[:15]),
        "mobility": {r["predicted_mobility"]: r["cnt"] for r in kpc_mob},
    }


def integron_ml_features():
    """
    Build feature matrix to predict integron carriage.
    Target: plasmid has qacEdelta1 + sul1 (class 1 integron).
    """
    # Positive: plasmids with integron
    positives = q("""
        SELECT DISTINCT a1.NUCCORE_ACC
        FROM amr a1
        JOIN amr a2 ON a1.NUCCORE_ACC = a2.NUCCORE_ACC
        WHERE a1.gene_symbol = 'qacEdelta1' AND a2.gene_symbol = 'sul1'
    """)
    pos_accs = {r["NUCCORE_ACC"] for r in positives}

    # Get features for all plasmids
    rows = q("""
        SELECT n.NUCCORE_ACC, n.NUCCORE_Length AS length, n.NUCCORE_GC AS gc,
               t.predicted_mobility AS mobility, t.rep_type,
               tx.TAXONOMY_genus AS genus
        FROM nuccore n
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        JOIN taxonomy tx ON n.TAXONOMY_UID = tx.TAXONOMY_UID
        WHERE t.predicted_mobility != '' AND tx.TAXONOMY_genus != ''
          AND n.NUCCORE_Length > 0
    """)

    for r in rows:
        r["has_integron"] = 1 if r["NUCCORE_ACC"] in pos_accs else 0

    return rows


def transposon_amr_cooccurrence():
    """Which IS families co-occur with which AMR drug classes."""
    import re as re_mod
    # Get IS family per plasmid from PGAP
    pgap = q("""
        SELECT NUCCORE_ACC, product FROM pgap_features
        WHERE category = 'TRANSPOSASE'
    """)
    plasmid_is = {}
    for r in pgap:
        m = re_mod.search(r"(IS\w+)\s+family", r["product"])
        if m:
            fam = m.group(1)
            plasmid_is.setdefault(r["NUCCORE_ACC"], set()).add(fam)

    # Get AMR drug class per plasmid
    amr_rows = q("""
        SELECT NUCCORE_ACC, drug_class FROM amr
        WHERE drug_class IS NOT NULL AND drug_class != ''
    """)
    plasmid_amr = {}
    for r in amr_rows:
        plasmid_amr.setdefault(r["NUCCORE_ACC"], set()).add(r["drug_class"])

    # Build co-occurrence matrix
    cooccur = {}
    for acc, is_fams in plasmid_is.items():
        amr_classes = plasmid_amr.get(acc, set())
        for is_fam in is_fams:
            for amr_cls in amr_classes:
                cooccur.setdefault(is_fam, {})
                cooccur[is_fam][amr_cls] = cooccur[is_fam].get(amr_cls, 0) + 1

    return cooccur


# ── Retro-mobilization / HGT routes ───────────────────────────────

def retromobilization_analysis():
    """
    Analyse transfer routes for non-mobilizable plasmids:
    1. Retro-mobilization (conjugative partner in same host → can pull DNA back)
    2. Co-mobilization relay (mobilizable partner)
    3. Conduction (co-integration via shared IS elements)
    4. Transduction (small plasmids fit in phage capsids)
    """
    total_nonmob = scalar(
        "SELECT COUNT(*) FROM typing WHERE predicted_mobility='non-mobilizable'"
    ) or 1

    # Retro-mobilization: non-mob + conjugative in same host
    retro = scalar("""
        SELECT COUNT(DISTINCT n2.NUCCORE_ACC)
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t1 ON n1.NUCCORE_ACC = t1.NUCCORE_ACC
        JOIN typing t2 ON n2.NUCCORE_ACC = t2.NUCCORE_ACC
        WHERE t1.predicted_mobility = 'conjugative'
          AND t2.predicted_mobility = 'non-mobilizable'
          AND n1.BIOSAMPLE_UID > 0
    """) or 0

    # Retro-mob with AMR
    retro_amr = scalar("""
        SELECT COUNT(DISTINCT n2.NUCCORE_ACC)
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t1 ON n1.NUCCORE_ACC = t1.NUCCORE_ACC
        JOIN typing t2 ON n2.NUCCORE_ACC = t2.NUCCORE_ACC
        JOIN amr a ON n2.NUCCORE_ACC = a.NUCCORE_ACC
        WHERE t1.predicted_mobility = 'conjugative'
          AND t2.predicted_mobility = 'non-mobilizable'
          AND n1.BIOSAMPLE_UID > 0
          AND a.drug_class IS NOT NULL AND a.drug_class != ''
    """) or 0

    # Mobilizable relay partner
    mob_relay = scalar("""
        SELECT COUNT(DISTINCT n2.NUCCORE_ACC)
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t1 ON n1.NUCCORE_ACC = t1.NUCCORE_ACC
        JOIN typing t2 ON n2.NUCCORE_ACC = t2.NUCCORE_ACC
        WHERE t1.predicted_mobility = 'mobilizable'
          AND t2.predicted_mobility = 'non-mobilizable'
          AND n1.BIOSAMPLE_UID > 0
    """) or 0

    # Small non-mob (<10kb) = transduction candidates
    small = scalar("""
        SELECT COUNT(*) FROM nuccore n
        JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE t.predicted_mobility = 'non-mobilizable' AND n.NUCCORE_Length < 10000
    """) or 0

    # Non-mob with IS elements (conduction potential via co-integration)
    conduction = scalar("""
        SELECT COUNT(DISTINCT p.NUCCORE_ACC) FROM pgap_features p
        JOIN typing t ON p.NUCCORE_ACC = t.NUCCORE_ACC
        WHERE t.predicted_mobility = 'non-mobilizable'
          AND p.category = 'TRANSPOSASE'
    """) or 0

    # Conjugative Inc groups in retro-mob pairs
    retro_inc = {}
    rows = q("""
        SELECT t1.rep_type, COUNT(DISTINCT n2.NUCCORE_ACC) as cnt
        FROM nuccore n1
        JOIN nuccore n2 ON n1.BIOSAMPLE_UID = n2.BIOSAMPLE_UID
                       AND n1.NUCCORE_ACC != n2.NUCCORE_ACC
        JOIN typing t1 ON n1.NUCCORE_ACC = t1.NUCCORE_ACC
        JOIN typing t2 ON n2.NUCCORE_ACC = t2.NUCCORE_ACC
        WHERE t1.predicted_mobility = 'conjugative'
          AND t2.predicted_mobility = 'non-mobilizable'
          AND n1.BIOSAMPLE_UID > 0 AND t1.rep_type != ''
        GROUP BY t1.rep_type ORDER BY cnt DESC LIMIT 50
    """)
    for r in rows:
        for rt in r["rep_type"].split(","):
            rt = rt.strip()
            if rt:
                retro_inc[rt] = retro_inc.get(rt, 0) + r["cnt"]

    return {
        "total_nonmob": total_nonmob,
        "retro_mob": retro,
        "retro_mob_amr": retro_amr,
        "mob_relay": mob_relay,
        "small_transduction": small,
        "conduction_potential": conduction,
        "retro_inc": dict(sorted(retro_inc.items(), key=lambda x: -x[1])[:12]),
        "routes": {
            "Retro-mobilization": retro,
            "Mobilizable relay": mob_relay,
            "Transduction (<10kb)": small,
            "Conduction (has IS/Tn)": conduction,
        },
    }


# ── Data export ────────────────────────────────────────────────────

EXPORT_QUERIES = {
    "all_plasmids": (
        "All Plasmids (PLSDB)",
        """SELECT n.NUCCORE_ACC, n.NUCCORE_Description, n.NUCCORE_Length,
                  n.NUCCORE_GC, n.NUCCORE_Topology, n.NUCCORE_CreateDate, n.NUCCORE_Source,
                  t.TAXONOMY_genus, t.TAXONOMY_species
           FROM nuccore n
           LEFT JOIN taxonomy t ON n.TAXONOMY_UID = t.TAXONOMY_UID""",
    ),
    "typing": (
        "MOB Typing & Inc Groups",
        """SELECT NUCCORE_ACC, rep_type, relaxase_type, mpf_type, orit_type,
                  predicted_mobility, predicted_host_range_overall_name,
                  PMLST_scheme, PMLST_alleles
           FROM typing""",
    ),
    "amr": (
        "AMR Annotations",
        """SELECT NUCCORE_ACC, gene_symbol, gene_name, drug_class,
                  antimicrobial_agent, input_gene_start, input_gene_stop,
                  strand_orientation, sequence_identity, coverage_percentage
           FROM amr""",
    ),
    "geography": (
        "Geographic Data",
        """SELECT n.NUCCORE_ACC, b.country, b.lat, b.lng,
                  b.host_category, b.ecosystem,
                  t.predicted_mobility, t.rep_type
           FROM nuccore n
           JOIN biosample_location b ON n.BIOSAMPLE_UID = b.BIOSAMPLE_UID
           JOIN typing t ON n.NUCCORE_ACC = t.NUCCORE_ACC
           WHERE b.country != ''""",
    ),
    "ncbi_extra": (
        "NCBI Extra Plasmids",
        "SELECT * FROM ncbi_extra",
    ),
}


def export_csv(dataset_key):
    """Export a dataset as CSV string."""
    import io, csv as csv_mod
    if dataset_key not in EXPORT_QUERIES:
        return None
    _, sql = EXPORT_QUERIES[dataset_key]
    rows = q(sql)
    if not rows:
        return None
    output = io.StringIO()
    writer = csv_mod.DictWriter(output, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


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
