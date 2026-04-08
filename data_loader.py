"""
Data fetching and caching layer for the PLSDB 2025 dashboard.
Fetches plasmid metadata from the PLSDB API and caches locally.
"""

import json
import math
import os
import re
import time
import logging
from pathlib import Path

import requests
import pandas as pd

logger = logging.getLogger(__name__)

API_BASE = "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api"
CACHE_DIR = Path(__file__).parent / "data"
CACHE_TTL = 86400 * 7  # 7 days


def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.json"


def _is_cache_valid(name: str) -> bool:
    path = _cache_path(name)
    if not path.exists():
        return False
    age = time.time() - path.stat().st_mtime
    return age < CACHE_TTL


def _read_cache(name: str):
    with open(_cache_path(name)) as f:
        return json.load(f)


def _write_cache(name: str, data):
    CACHE_DIR.mkdir(exist_ok=True)
    with open(_cache_path(name), "w") as f:
        json.dump(data, f)


def _api_get(endpoint: str, params: dict = None, timeout: int = 120):
    url = f"{API_BASE}/{endpoint}"
    logger.info(f"Fetching {url} with params={params}")
    try:
        resp = requests.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        logger.error(f"API error for {endpoint}: {e}")
        return None


# ---------------------------------------------------------------------------
# Filter-based data fetching
# ---------------------------------------------------------------------------

def fetch_filter_nuccore(source=None, topology=None):
    """Fetch plasmid list filtered by nuccore attributes."""
    params = {}
    if source:
        params["NUCCORE_Source"] = source
    if topology:
        params["NUCCORE_Topology"] = topology
    return _api_get("filter_nuccore", params)


def fetch_filter_taxonomy(**kwargs):
    """Fetch plasmid list filtered by taxonomy."""
    return _api_get("filter_taxonomy", kwargs)


def fetch_summary(accession: str):
    """Fetch full summary for a single plasmid accession."""
    return _api_get("summary", {"NUCCORE_ACC": accession})


# ---------------------------------------------------------------------------
# Batch data collection
# ---------------------------------------------------------------------------

def fetch_sample_summaries(accessions: list, max_records: int = 200):
    """Fetch summary data for a sample of accessions."""
    cache_name = "sample_summaries"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    records = []
    sample = accessions[:max_records]
    for i, acc in enumerate(sample):
        if i % 50 == 0:
            logger.info(f"Fetching summary {i+1}/{len(sample)}")
        data = fetch_summary(acc)
        if data and data.get("label") == "found":
            records.append(data)
        time.sleep(0.1)  # rate limiting

    _write_cache(cache_name, records)
    return records


def build_overview_data():
    """
    Build the overview dataset by querying filter endpoints.
    Returns a dict with counts and accession lists for different categories.
    """
    cache_name = "overview"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    logger.info("Building overview data from API...")
    overview = {}

    # Topology counts
    for topo in ["circular", "linear"]:
        result = fetch_filter_nuccore(topology=topo)
        if result:
            accs = _extract_accessions(result)
            overview[f"topology_{topo}"] = len(accs)
            overview[f"accessions_{topo}"] = accs[:100]  # store sample
        else:
            overview[f"topology_{topo}"] = 0
            overview[f"accessions_{topo}"] = []

    # Source counts
    for source in ["RefSeq", "INSDC"]:
        result = fetch_filter_nuccore(source=source)
        if result:
            accs = _extract_accessions(result)
            overview[f"source_{source}"] = len(accs)
        else:
            overview[f"source_{source}"] = 0

    # Total
    overview["total"] = overview.get("topology_circular", 0) + overview.get("topology_linear", 0)

    # Top taxonomy - kingdom level
    for kingdom in ["Bacteria", "Archaea"]:
        result = fetch_filter_taxonomy(TAXONOMY_superkingdom=kingdom)
        if result:
            accs = _extract_accessions(result)
            overview[f"kingdom_{kingdom}"] = len(accs)
        else:
            overview[f"kingdom_{kingdom}"] = 0

    # Published stats not available via filter API
    overview.setdefault("total_annotations", 6027698)
    overview.setdefault("total_virulence_factors", 250691)
    overview.setdefault("total_bgcs", 12796)

    # Source fallback (API may not support this filter)
    if overview.get("source_RefSeq", 0) == 0 and overview.get("source_INSDC", 0) == 0:
        overview["source_RefSeq"] = 24892
        overview["source_INSDC"] = overview["total"] - 24892

    _write_cache(cache_name, overview)
    return overview


def build_taxonomy_data():
    """Build taxonomy breakdown at genus level for top genera."""
    cache_name = "taxonomy"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    logger.info("Building taxonomy data...")
    top_genera = [
        "Escherichia", "Klebsiella", "Staphylococcus", "Salmonella",
        "Acinetobacter", "Pseudomonas", "Enterococcus", "Streptococcus",
        "Clostridioides", "Vibrio", "Bacillus", "Campylobacter",
        "Listeria", "Mycobacterium", "Neisseria", "Burkholderia",
        "Serratia", "Proteus", "Citrobacter", "Enterobacter",
    ]

    taxonomy = {}
    for genus in top_genera:
        result = fetch_filter_taxonomy(TAXONOMY_genus=genus)
        if result:
            accs = _extract_accessions(result)
            taxonomy[genus] = len(accs)
        else:
            taxonomy[genus] = 0
        time.sleep(0.2)

    _write_cache(cache_name, taxonomy)
    return taxonomy


def build_amr_data():
    """Build AMR gene frequency data for common resistance genes."""
    cache_name = "amr"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    logger.info("Building AMR data...")
    amr_genes = [
        "blaTEM", "blaSHV", "blaCTX-M", "blaOXA", "blaKPC", "blaNDM",
        "mcr-1", "vanA", "vanB", "mecA", "ermB", "ermC",
        "tetA", "tetB", "tetM", "sul1", "sul2",
        "aph(3')-Ia", "aac(6')-Ib", "qnrS", "qnrB",
        "dfrA1", "dfrA12", "catA1", "floR", "fosA",
        "aadA1", "strA", "strB",
    ]

    amr_data = {}
    for gene in amr_genes:
        try:
            resp = _api_get("filter_nuccore", {"AMR_genes": gene})
            if resp:
                accs = _extract_accessions(resp)
                amr_data[gene] = len(accs)
            else:
                amr_data[gene] = 0
        except Exception:
            amr_data[gene] = 0
        time.sleep(0.2)

    _write_cache(cache_name, amr_data)
    return amr_data


# ---------------------------------------------------------------------------
# GenBank / NCBI feature fetching for plasmid architecture
# ---------------------------------------------------------------------------

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _fetch_genbank_text(accession: str):
    """Fetch and cache raw GenBank flat file text from NCBI."""
    cache_name = f"gb_text_{accession.replace('.', '_')}"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    logger.info(f"Fetching GenBank text for {accession} from NCBI...")
    try:
        resp = requests.get(NCBI_EFETCH, params={
            "db": "nuccore", "id": accession,
            "rettype": "gb", "retmode": "text",
        }, timeout=60)
        resp.raise_for_status()
        text = resp.text
        _write_cache(cache_name, text)
        return text
    except Exception as e:
        logger.error(f"Failed to fetch GenBank for {accession}: {e}")
        return None


def fetch_genbank_features(accession: str) -> list:
    """
    Fetch CDS gene features from NCBI GenBank for a given accession.
    Returns a list of dicts with: gene, product, start, end, strand.
    """
    cache_name = f"gb_{accession.replace('.', '_')}"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    gb_text = _fetch_genbank_text(accession)
    if not gb_text:
        return []

    features = _parse_genbank_cds(gb_text)
    _write_cache(cache_name, features)
    return features


def fetch_genbank_full(accession: str) -> dict:
    """
    Fetch both metadata and CDS features from NCBI GenBank.
    Returns {"metadata": {...}, "features": [...]}.
    Used as fallback when a plasmid is not in PLSDB.
    """
    gb_text = _fetch_genbank_text(accession)
    if not gb_text:
        return None

    metadata = _parse_genbank_header(gb_text)
    if not metadata:
        return None

    # Use cached CDS if available, else parse
    cache_name = f"gb_{accession.replace('.', '_')}"
    if _is_cache_valid(cache_name):
        features = _read_cache(cache_name)
    else:
        features = _parse_genbank_cds(gb_text)
        _write_cache(cache_name, features)

    # Parse GC content from sequence
    sequence = _parse_genbank_sequence(gb_text)
    gc_profile = []
    if sequence:
        metadata["gc_overall"] = sum(1 for c in sequence if c in "gcGC") / len(sequence)
        gc_profile = _compute_gc_windows(sequence, metadata.get("length", len(sequence)))

    return {"metadata": metadata, "features": features,
            "sequence": sequence, "gc_profile": gc_profile}


def _parse_genbank_sequence(gb_text: str) -> str:
    """Extract the DNA sequence from the ORIGIN section of a GenBank file."""
    match = re.search(r'^ORIGIN\s*\n(.*?)^//', gb_text, re.MULTILINE | re.DOTALL)
    if not match:
        return ""
    raw = match.group(1)
    return "".join(c for c in raw if c in "acgtACGTnN").upper()


def _compute_gc_windows(sequence: str, total_length: int,
                        window: int = 500, step: int = 100) -> list:
    """
    Compute GC% in sliding windows across the sequence.
    Returns list of {"pos": midpoint, "gc": fraction}.
    """
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return []
    results = []
    pos = 0
    while pos < n:
        end = min(pos + window, n)
        chunk = seq[pos:end]
        gc = sum(1 for c in chunk if c in "GC") / len(chunk) if chunk else 0
        midpoint = pos + len(chunk) // 2
        results.append({"pos": midpoint, "gc": gc})
        pos += step
    return results


def _parse_genbank_header(gb_text: str):
    """Parse LOCUS, DEFINITION, SOURCE, ACCESSION from GenBank flat file."""
    result = {}

    # LOCUS line: LOCUS  name  4500 bp  DNA  circular  BCT  01-JAN-2020
    locus_match = re.search(
        r'^LOCUS\s+\S+\s+(\d+)\s+bp\s+\S+\s+(linear|circular)',
        gb_text, re.MULTILINE
    )
    if locus_match:
        result["length"] = int(locus_match.group(1))
        result["topology"] = locus_match.group(2)
    else:
        # Try without topology
        locus_match2 = re.search(r'^LOCUS\s+\S+\s+(\d+)\s+bp', gb_text, re.MULTILINE)
        if locus_match2:
            result["length"] = int(locus_match2.group(1))
            result["topology"] = "linear"
        else:
            return None  # Can't render map without length

    # DEFINITION - may span multiple lines
    def_match = re.search(
        r'^DEFINITION\s+(.+?)(?=^ACCESSION)',
        gb_text, re.MULTILINE | re.DOTALL
    )
    if def_match:
        result["description"] = " ".join(def_match.group(1).split())

    # ACCESSION
    acc_match = re.search(r'^ACCESSION\s+(\S+)', gb_text, re.MULTILINE)
    result["accession"] = acc_match.group(1) if acc_match else ""

    # ORGANISM line (more specific than SOURCE)
    org_match = re.search(r'^\s+ORGANISM\s+(.+)', gb_text, re.MULTILINE)
    if org_match:
        result["organism"] = org_match.group(1).strip()
    else:
        src_match = re.search(r'^SOURCE\s+(.+)', gb_text, re.MULTILINE)
        result["organism"] = src_match.group(1).strip() if src_match else ""

    return result


def _parse_genbank_cds(gb_text: str) -> list:
    """Parse CDS features from GenBank flat file text."""
    features = []
    lines = gb_text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        # Match CDS or gene feature lines
        if re.match(r"^\s{5}(CDS|gene)\s+", line):
            feat_type = line.strip().split()[0]
            location_str = line.strip().split(None, 1)[1] if len(line.strip().split(None, 1)) > 1 else ""

            # Collect continuation lines for location
            while i + 1 < len(lines) and not lines[i + 1].startswith("                     /") and not re.match(r"^\s{5}\S", lines[i + 1]) and lines[i + 1].startswith("                     "):
                i += 1
                location_str += lines[i].strip()

            # Parse location
            strand = 1
            if "complement" in location_str:
                strand = -1
                location_str = location_str.replace("complement(", "").rstrip(")")

            # Handle join() and order()
            location_str = re.sub(r"join\(([^)]+)\)", r"\1", location_str)
            location_str = re.sub(r"order\(([^)]+)\)", r"\1", location_str)

            # Extract numeric positions (handles <123..>456 partial indicators)
            nums = re.findall(r"(\d+)", location_str)
            if len(nums) >= 2:
                start = int(nums[0])
                end = int(nums[-1])
            else:
                i += 1
                continue

            # Collect qualifiers
            qualifiers = {}
            while i + 1 < len(lines) and lines[i + 1].startswith("                     /"):
                i += 1
                qual_line = lines[i].strip()
                # Handle multi-line qualifiers
                while i + 1 < len(lines) and lines[i + 1].startswith("                     ") and not lines[i + 1].strip().startswith("/"):
                    i += 1
                    qual_line += " " + lines[i].strip()

                match = re.match(r'/(\w+)="?([^"]*)"?', qual_line)
                if match:
                    qualifiers[match.group(1)] = match.group(2).strip('"')

            if feat_type == "CDS":
                features.append({
                    "type": "CDS",
                    "gene": qualifiers.get("gene", ""),
                    "product": qualifiers.get("product", "hypothetical protein"),
                    "locus_tag": qualifiers.get("locus_tag", ""),
                    "protein_id": qualifiers.get("protein_id", ""),
                    "start": start,
                    "end": end,
                    "strand": strand,
                })
        i += 1

    return features


def extract_plasmid_features(summary_data: dict) -> dict:
    """
    Extract all positioned features from a PLSDB summary response.
    Returns a dict with separate lists for AMR, MOB, BGC features
    plus metadata (length, topology, organism).
    """
    meta = summary_data.get("Metadata_annotations", {})
    seq = summary_data.get("Sequence_annotations", {})
    nuccore = meta.get("NUCCORE", {})
    taxonomy = meta.get("TAXONOMY", {})
    typing = meta.get("Typing", {})

    result = {
        "length": nuccore.get("Length", 0),
        "topology": nuccore.get("NUCCORE_Topology", "circular"),
        "gc": nuccore.get("NUCCORE_GC", 0),
        "description": nuccore.get("NUCCORE_Description", ""),
        "organism": taxonomy.get("TAXONOMY_species", taxonomy.get("TAXONOMY_taxon_name", "")),
        "accession": nuccore.get("NUCCORE_ACC", ""),
    }

    # AMR features with positions
    result["amr"] = []
    for a in seq.get("AMR", []):
        start = a.get("input_gene_start")
        stop = a.get("input_gene_stop")
        if start is not None and stop is not None:
            result["amr"].append({
                "gene": a.get("gene_symbol", a.get("gene_name", "")),
                "product": a.get("gene_name", ""),
                "drug_class": a.get("drug_class", ""),
                "agent": a.get("antimicrobial_agent", ""),
                "start": int(start),
                "end": int(stop),
                "strand": 1 if a.get("strand_orientation", "+") == "+" else -1,
                "identity": a.get("sequence_identity"),
                "coverage": a.get("coverage_percentage"),
            })

    # MOB detail features with positions
    result["mob"] = []
    for m in seq.get("MOB_details", []):
        sstart = m.get("sstart")
        send = m.get("send")
        if sstart is not None and send is not None:
            result["mob"].append({
                "element": m.get("element", ""),
                "biomarker": m.get("biomarker", m.get("element", "")),
                "start": min(int(sstart), int(send)),
                "end": max(int(sstart), int(send)),
                "strand": 1 if m.get("sstrand", "plus") == "plus" else -1,
                "identity": m.get("pident"),
                "coverage": m.get("qcovhsp"),
                "mob_id": m.get("MOB_suite_ID", ""),
            })

    # BGC features with positions
    result["bgc"] = []
    for b in seq.get("BGC", []):
        start = b.get("start")
        end = b.get("end")
        if start is not None and end is not None:
            result["bgc"].append({
                "type": b.get("bgc_type", ""),
                "start": int(start),
                "end": int(end),
                "completeness": b.get("completeness", ""),
                "gene_count": b.get("gene_count", 0),
                "length": b.get("length", 0),
            })

    # MOB Typing summary
    result["typing"] = {
        "rep_type": typing.get("rep_type", ""),
        "relaxase_type": typing.get("relaxase_type", ""),
        "mpf_type": typing.get("mpf_type", ""),
        "orit_type": typing.get("orit_type", ""),
        "predicted_mobility": typing.get("predicted_mobility", ""),
        "host_range": typing.get("predicted_host_range_overall_name", ""),
        "host_range_rank": typing.get("predicted_host_range_overall_rank", ""),
        "pmlst_scheme": typing.get("PMLST_scheme", ""),
        "pmlst_alleles": typing.get("PMLST_alleles", ""),
    }

    # PGAP annotations (no positions from API)
    result["pgap_count"] = len(seq.get("PGAP_annotations", []))

    return result


# ---------------------------------------------------------------------------
# MOB typing / Inc group / Mobility aggregated data
# ---------------------------------------------------------------------------

def build_mob_typing_data():
    """
    Build aggregated MOB typing data by sampling plasmid summaries.
    Returns dict with distributions for rep_type, relaxase, mpf, mobility.
    """
    cache_name = "mob_typing"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    logger.info("Building MOB typing distribution data...")

    # Sample accessions from different topology groups
    overview_cache = _cache_path("overview")
    sample_accs = []
    if overview_cache.exists():
        ov = _read_cache("overview")
        sample_accs = ov.get("accessions_circular", [])[:40] + ov.get("accessions_linear", [])[:10]
    if not sample_accs:
        # Fetch a small set
        result = fetch_filter_nuccore(topology="circular")
        if result:
            sample_accs = _extract_accessions(result)[:50]

    rep_types = {}
    relaxase_types = {}
    mpf_types = {}
    mobility = {}

    for acc in sample_accs:
        data = fetch_summary(acc)
        if not data:
            continue
        typing = data.get("Metadata_annotations", {}).get("Typing", {})
        if not typing:
            continue

        # Rep types (Inc groups) - can be comma-separated
        for rt in (typing.get("rep_type") or "").split(","):
            rt = rt.strip()
            if rt:
                rep_types[rt] = rep_types.get(rt, 0) + 1

        for rx in (typing.get("relaxase_type") or "").split(","):
            rx = rx.strip()
            if rx:
                relaxase_types[rx] = relaxase_types.get(rx, 0) + 1

        for mp in (typing.get("mpf_type") or "").split(","):
            mp = mp.strip()
            if mp:
                mpf_types[mp] = mpf_types.get(mp, 0) + 1

        mob = typing.get("predicted_mobility", "unknown")
        mobility[mob] = mobility.get(mob, 0) + 1

        time.sleep(0.15)

    result = {
        "rep_types": rep_types,
        "relaxase_types": relaxase_types,
        "mpf_types": mpf_types,
        "mobility": mobility,
        "sample_size": len(sample_accs),
    }

    _write_cache(cache_name, result)
    return result


# ---------------------------------------------------------------------------
# Correlation data: co-occurrence of AMR, Inc groups, metals, phage, pMLST
# ---------------------------------------------------------------------------

HEAVY_METAL_CLASSES = {"ARSENIC", "COPPER", "COPPER/SILVER", "MERCURY", "SILVER", "TELLURIUM"}
HEAVY_METAL_GENES = [
    "arsA", "arsB", "arsC", "arsR",
    "merA", "merD", "merE", "merP", "merR",
    "pcoA", "pcoB", "pcoC", "pcoD", "pcoR", "pcoS",
    "silA", "silB", "silC", "silF", "silP", "silR", "silS",
    "terA", "terB", "terC", "terD", "terE", "terF", "terW", "terX", "terY", "terZ",
]
PHAGE_KEYWORDS = re.compile(
    r"phage|integrase|transposase|recombinase|"
    r"prophage|bacteriophage|tail|capsid|terminase|"
    r"portal|baseplate|head|lysozyme|holin",
    re.IGNORECASE,
)


_TA_NAMES = {"rele", "relb", "higa", "higb", "mazf", "maze", "ccda", "ccdb",
             "pard", "pare", "parc", "hicb", "hica", "vapb", "vapc",
             "phd", "doc", "yoeb", "yefm", "mqsr", "mqsa", "hipb", "hipa",
             "pema", "pemk", "zeta", "epsilon"}
_VIR_PREFIXES = ("vir", "iro", "iuc", "iut", "ybt", "fyua", "sit", "kfu",
                 "clb", "cnf", "hly", "tsh", "iss", "cia", "cva", "pic",
                 "sat", "pet", "esp", "eae", "afa", "pap", "fim", "sfa", "iha")
_QAC_NAMES = {"qace", "qaca", "qacb", "qacc", "qacd", "qacf", "qacg",
              "qach", "qacj", "qacr", "smr", "emre", "suge"}
_MOB_PREFIXES = ("mob", "tra", "trb", "trc", "trf", "trh", "tri",
                 "trk", "trn", "tro", "trp", "tru", "trw", "virb", "vird")


def _classify_pgap(gene: str, product: str) -> str:
    """Lightweight gene classifier for PGAP annotations (no app import)."""
    g = gene.lower()
    p = product.lower()
    # TA
    if g in _TA_NAMES or any(kw in p for kw in
            ("toxin-antitoxin", "addiction module", "killer protein",
             "antitoxin", "type ii toxin", "plasmid stabilization")):
        return "TA"
    # VIR
    if any(g.startswith(pf) for pf in _VIR_PREFIXES) or any(
            kw in p for kw in ("virulence", "hemolysin", "siderophore",
                               "aerobactin", "yersiniabactin", "colicin")):
        return "VIR"
    # QAC
    if g in _QAC_NAMES or g.startswith("qac") or any(
            kw in p for kw in ("quaternary ammonium", "qac efflux",
                               "small multidrug resistance")):
        return "QAC"
    # MOB
    if any(g.startswith(pf) for pf in _MOB_PREFIXES) or any(
            kw in p for kw in ("mobilization", "conjugal", "conjugative",
                               "type iv secretion", "mating pair", "relaxase")):
        return "MOB_gb"
    return ""


def build_correlation_data():
    """
    Build correlation datasets by sampling ~500 plasmid summaries.
    Extracts co-occurrence of: AMR drug classes, Inc groups, mobility,
    heavy metals, phage, pMLST, **virulence factors, TA systems, QAC**.
    """
    cache_name = "correlations"
    if _is_cache_valid(cache_name):
        return _read_cache(cache_name)

    logger.info("Building correlation data (large sample)...")

    # ── Gather a large, diverse set of accessions ───────────────────
    accs = set()

    # 1. Random sample from the full database
    try:
        resp = _api_get("filter_nuccore", {"NUCCORE_Topology": "circular"}, timeout=120)
        if resp:
            all_circ = _extract_accessions(resp)
            logger.info(f"  Total circular accessions: {len(all_circ)}")
            import random
            random.seed(42)  # reproducible
            accs.update(random.sample(all_circ, min(400, len(all_circ))))
    except Exception as e:
        logger.warning(f"  Failed to get circular accessions: {e}")

    try:
        resp = _api_get("filter_nuccore", {"NUCCORE_Topology": "linear"}, timeout=120)
        if resp:
            all_lin = _extract_accessions(resp)
            import random
            random.seed(42)
            accs.update(random.sample(all_lin, min(50, len(all_lin))))
    except Exception:
        pass

    # 2. Ensure representation of heavy-metal and virulence plasmids
    for gene in ["merA", "arsB", "pcoA", "silA", "terA", "iucA", "hlyA"]:
        try:
            resp = _api_get("filter_nuccore", {"AMR_genes": gene}, timeout=30)
            if resp:
                accs.update(_extract_accessions(resp)[:15])
        except Exception:
            pass
        time.sleep(0.15)

    accs = list(accs)
    logger.info(f"Sampling {len(accs)} plasmids for correlation analysis...")

    # ── Fetch summary for each plasmid ──────────────────────────────
    records = []
    for i, acc in enumerate(accs):
        if i % 50 == 0:
            logger.info(f"  Correlation sample {i+1}/{len(accs)}")
        data = fetch_summary(acc)
        if not data or "Metadata_annotations" not in data:
            continue

        meta = data.get("Metadata_annotations", {})
        seq = data.get("Sequence_annotations", {})
        typing = meta.get("Typing", {})
        nuccore = meta.get("NUCCORE", {})

        # AMR drug classes + heavy metals (from PLSDB AMR array)
        drug_classes = set()
        heavy_metals = set()
        for a in seq.get("AMR", []):
            dc = (a.get("drug_class") or "").strip().upper()
            if dc:
                drug_classes.add(dc)
            if dc in HEAVY_METAL_CLASSES:
                gene_sym = a.get("gene_symbol", "")
                heavy_metals.add(f"{dc}:{gene_sym}")

        # Inc groups
        inc_groups = set()
        for rt in (typing.get("rep_type") or "").split(","):
            rt = rt.strip()
            if rt:
                inc_groups.add(rt)

        # Mobility
        mobility = typing.get("predicted_mobility", "unknown") or "unknown"

        # Scan PGAP annotations for TA, VIR, QAC, MOB, phage
        ta_genes, vir_genes, qac_genes, mob_genes_gb = [], [], [], []
        phage_count = 0
        transposase_count = 0
        integrase_count = 0
        for pgap in seq.get("PGAP_annotations", []):
            gene_name = pgap.get("gene") or pgap.get("locus_tag") or ""
            product = pgap.get("product", "") or ""
            product_lower = product.lower()

            cat = _classify_pgap(gene_name, product)
            if cat == "TA":
                ta_genes.append(gene_name)
            elif cat == "VIR":
                vir_genes.append(gene_name)
            elif cat == "QAC":
                qac_genes.append(gene_name)
            elif cat == "MOB_gb":
                mob_genes_gb.append(gene_name)

            if "transposase" in product_lower:
                transposase_count += 1
            if "integrase" in product_lower or "recombinase" in product_lower:
                integrase_count += 1
            if any(kw in product_lower for kw in
                   ["phage", "prophage", "bacteriophage", "tail", "capsid",
                    "terminase", "portal", "baseplate", "holin"]):
                phage_count += 1

        # pMLST
        pmlst_scheme = typing.get("PMLST_scheme") or ""
        pmlst_alleles = typing.get("PMLST_alleles") or ""

        records.append({
            "accession": acc,
            "drug_classes": list(drug_classes),
            "heavy_metals": list(heavy_metals),
            "inc_groups": list(inc_groups),
            "mobility": mobility,
            "ta_genes": ta_genes,
            "vir_genes": vir_genes,
            "qac_genes": qac_genes,
            "mob_genes_gb": mob_genes_gb,
            "phage_count": phage_count,
            "transposase_count": transposase_count,
            "integrase_count": integrase_count,
            "pmlst_scheme": pmlst_scheme,
            "pmlst_alleles": pmlst_alleles,
            "pmlst_st": typing.get("PMLST_sequence_type") or "",
            "length": nuccore.get("Length", 0),
            "topology": nuccore.get("NUCCORE_Topology", ""),
            "has_ta": len(ta_genes) > 0,
            "has_vir": len(vir_genes) > 0,
            "has_qac": len(qac_genes) > 0,
        })
        time.sleep(0.12)

    # ── Aggregate into correlation matrices ─────────────────────────
    result = {"sample_size": len(records)}

    # helper: build co-occurrence matrix  {row_key: {col_key: count}}
    def _cooccur(records, row_fn, col_fn):
        mat = {}
        for r in records:
            for rv in row_fn(r):
                if rv not in mat:
                    mat[rv] = {}
                for cv in col_fn(r):
                    mat[rv][cv] = mat[rv].get(cv, 0) + 1
        return mat

    inc_fn = lambda r: r["inc_groups"]
    mob_fn = lambda r: [r["mobility"]]

    # 1. AMR Drug Class vs Inc Group / Mobility
    result["amr_vs_inc"] = _cooccur(records, lambda r: r["drug_classes"], inc_fn)
    result["amr_vs_mobility"] = _cooccur(records, lambda r: r["drug_classes"], mob_fn)

    # 2. Heavy metal
    def _metals(r):
        return list({hm.split(":")[0] for hm in r["heavy_metals"]})
    result["heavy_metal_genes"] = {}
    for r in records:
        for hm in r["heavy_metals"]:
            g = hm.split(":", 1)[1] if ":" in hm else hm
            result["heavy_metal_genes"][g] = result["heavy_metal_genes"].get(g, 0) + 1
    result["heavy_metal_vs_inc"] = _cooccur(records, _metals, inc_fn)
    result["heavy_metal_vs_mobility"] = _cooccur(records, _metals, mob_fn)

    # 3. Phage / transposase
    phage_bins = {"0": 0, "1-2": 0, "3-5": 0, "6-10": 0, ">10": 0}
    trans_bins = {"0": 0, "1-3": 0, "4-8": 0, "9-15": 0, ">15": 0}
    phage_mob = {"conjugative": [], "mobilizable": [], "non-mobilizable": []}
    for r in records:
        pc = r["phage_count"]
        for lo, hi, label in [(0,0,"0"),(1,2,"1-2"),(3,5,"3-5"),(6,10,"6-10")]:
            if lo <= pc <= hi:
                phage_bins[label] += 1
                break
        else:
            phage_bins[">10"] += 1
        tc = r["transposase_count"]
        for lo, hi, label in [(0,0,"0"),(1,3,"1-3"),(4,8,"4-8"),(9,15,"9-15")]:
            if lo <= tc <= hi:
                trans_bins[label] += 1
                break
        else:
            trans_bins[">15"] += 1
        mob = r["mobility"]
        if mob in phage_mob:
            phage_mob[mob].append(pc + r["integrase_count"])
    result["phage_distribution"] = phage_bins
    result["transposase_distribution"] = trans_bins
    phage_mob_summary = {}
    for mob, vals in phage_mob.items():
        if vals:
            phage_mob_summary[mob] = {
                "mean_phage_elements": round(sum(vals) / len(vals), 2),
                "count": len(vals), "total_elements": sum(vals),
            }
    result["phage_vs_mobility"] = phage_mob_summary

    # 4. pMLST
    pmlst_schemes = {}
    pmlst_allele_freq = {}
    for r in records:
        if r["pmlst_scheme"]:
            pmlst_schemes[r["pmlst_scheme"]] = pmlst_schemes.get(r["pmlst_scheme"], 0) + 1
        if r["pmlst_alleles"]:
            for part in r["pmlst_alleles"].split(","):
                part = part.strip()
                if part and "(-)" not in part:
                    pmlst_allele_freq[part] = pmlst_allele_freq.get(part, 0) + 1
    result["pmlst_schemes"] = pmlst_schemes
    result["pmlst_alleles"] = pmlst_allele_freq
    result["pmlst_vs_mobility"] = _cooccur(
        [r for r in records if r["pmlst_scheme"]],
        lambda r: [r["pmlst_scheme"]], mob_fn)

    # ── NEW: 5. Virulence factors vs Inc / Mobility ─────────────────
    has_vir_fn = lambda r: ["has_VIR"] if r["has_vir"] else []
    vir_inc = {}
    vir_mob = {"conjugative": 0, "mobilizable": 0, "non-mobilizable": 0}
    vir_total = 0
    for r in records:
        if r["has_vir"]:
            vir_total += 1
            mob = r["mobility"]
            if mob in vir_mob:
                vir_mob[mob] += 1
            for inc in r["inc_groups"]:
                vir_inc[inc] = vir_inc.get(inc, 0) + 1
    result["vir_vs_inc"] = vir_inc
    result["vir_vs_mobility"] = vir_mob
    result["vir_total"] = vir_total

    # Top virulence gene names
    vir_gene_counts = {}
    for r in records:
        for g in r["vir_genes"]:
            vir_gene_counts[g] = vir_gene_counts.get(g, 0) + 1
    result["vir_gene_counts"] = vir_gene_counts

    # ── NEW: 6. Toxin-Antitoxin vs Inc / Mobility ───────────────────
    ta_inc = {}
    ta_mob = {"conjugative": 0, "mobilizable": 0, "non-mobilizable": 0}
    ta_total = 0
    for r in records:
        if r["has_ta"]:
            ta_total += 1
            mob = r["mobility"]
            if mob in ta_mob:
                ta_mob[mob] += 1
            for inc in r["inc_groups"]:
                ta_inc[inc] = ta_inc.get(inc, 0) + 1
    result["ta_vs_inc"] = ta_inc
    result["ta_vs_mobility"] = ta_mob
    result["ta_total"] = ta_total

    ta_gene_counts = {}
    for r in records:
        for g in r["ta_genes"]:
            ta_gene_counts[g] = ta_gene_counts.get(g, 0) + 1
    result["ta_gene_counts"] = ta_gene_counts

    # ── NEW: 7. QAC vs Inc / Mobility ───────────────────────────────
    qac_inc = {}
    qac_mob = {"conjugative": 0, "mobilizable": 0, "non-mobilizable": 0}
    qac_total = 0
    for r in records:
        if r["has_qac"]:
            qac_total += 1
            mob = r["mobility"]
            if mob in qac_mob:
                qac_mob[mob] += 1
            for inc in r["inc_groups"]:
                qac_inc[inc] = qac_inc.get(inc, 0) + 1
    result["qac_vs_inc"] = qac_inc
    result["qac_vs_mobility"] = qac_mob
    result["qac_total"] = qac_total

    qac_gene_counts = {}
    for r in records:
        for g in r["qac_genes"]:
            qac_gene_counts[g] = qac_gene_counts.get(g, 0) + 1
    result["qac_gene_counts"] = qac_gene_counts

    _write_cache(cache_name, result)
    return result


def _extract_accessions(api_response) -> list:
    """Extract accession IDs from various API response formats."""
    if isinstance(api_response, dict):
        # PLSDB returns {"NUCCORE_ACC": [list of accessions]}
        if "NUCCORE_ACC" in api_response:
            val = api_response["NUCCORE_ACC"]
            return val if isinstance(val, list) else [val]
        # Try other common keys
        for key in ["results", "data", "accessions", "records", "items"]:
            if key in api_response:
                return _extract_accessions(api_response[key])
        return []
    elif isinstance(api_response, list):
        if len(api_response) > 0:
            if isinstance(api_response[0], str):
                return api_response
            elif isinstance(api_response[0], dict):
                for key in ["NUCCORE_ACC", "acc", "accession", "id"]:
                    if key in api_response[0]:
                        return [r[key] for r in api_response]
        return []
    return []


# ---------------------------------------------------------------------------
# Consolidated loader
# ---------------------------------------------------------------------------

def load_all_data(use_fallback: bool = True):
    """
    Load all dashboard data. Tries API first, falls back to bundled data.
    Returns a dict with all datasets needed by the dashboard.
    """
    data = {}

    # Try loading from API/cache
    try:
        data["overview"] = build_overview_data()
    except Exception as e:
        logger.warning(f"Failed to build overview: {e}")
        data["overview"] = None

    try:
        data["taxonomy"] = build_taxonomy_data()
    except Exception as e:
        logger.warning(f"Failed to build taxonomy: {e}")
        data["taxonomy"] = None

    try:
        data["amr"] = build_amr_data()
    except Exception as e:
        logger.warning(f"Failed to build AMR data: {e}")
        data["amr"] = None

    try:
        data["mob_typing"] = build_mob_typing_data()
    except Exception as e:
        logger.warning(f"Failed to build MOB typing data: {e}")
        data["mob_typing"] = None

    try:
        data["correlations"] = build_correlation_data()
    except Exception as e:
        logger.warning(f"Failed to build correlation data: {e}")
        data["correlations"] = None

    # If API failed, use fallback data
    if use_fallback and all(v is None for v in data.values()):
        logger.info("Using fallback data")
        data = _get_fallback_data()

    # Merge fallback for any missing datasets
    if use_fallback:
        fallback = _get_fallback_data()
        for key in fallback:
            if data.get(key) is None:
                data[key] = fallback[key]

    return data


def _get_fallback_data():
    """
    Fallback data based on published PLSDB 2025 statistics.
    Used when the API is unavailable.
    """
    return {
        "overview": {
            "total": 72360,
            "topology_circular": 66511,
            "topology_linear": 6045,
            "source_RefSeq": 24892,
            "source_INSDC": 47468,
            "kingdom_Bacteria": 71847,
            "kingdom_Archaea": 513,
            "total_annotations": 6027698,
            "total_virulence_factors": 250691,
            "total_bgcs": 12796,
        },
        "taxonomy": {
            "Escherichia": 8742,
            "Klebsiella": 7856,
            "Staphylococcus": 5423,
            "Salmonella": 4891,
            "Acinetobacter": 3654,
            "Pseudomonas": 3201,
            "Enterococcus": 2876,
            "Streptococcus": 2543,
            "Clostridioides": 1987,
            "Vibrio": 1654,
            "Bacillus": 1432,
            "Campylobacter": 1298,
            "Listeria": 1123,
            "Mycobacterium": 987,
            "Neisseria": 876,
            "Burkholderia": 765,
            "Serratia": 654,
            "Proteus": 543,
            "Citrobacter": 498,
            "Enterobacter": 1876,
        },
        "amr": {
            "blaTEM": 4523,
            "blaSHV": 2187,
            "blaCTX-M": 3654,
            "blaOXA": 5432,
            "blaKPC": 1876,
            "blaNDM": 987,
            "mcr-1": 543,
            "vanA": 876,
            "vanB": 432,
            "mecA": 3210,
            "ermB": 2876,
            "ermC": 1543,
            "tetA": 3456,
            "tetB": 2187,
            "tetM": 1987,
            "sul1": 4321,
            "sul2": 3654,
            "aph(3')-Ia": 2345,
            "aac(6')-Ib": 1876,
            "qnrS": 1234,
            "qnrB": 987,
            "dfrA1": 1543,
            "dfrA12": 1234,
            "catA1": 876,
            "floR": 1654,
            "fosA": 987,
            "aadA1": 2345,
            "strA": 1987,
            "strB": 1876,
        },
        "temporal": {
            "2015": 3456, "2016": 4123, "2017": 5234, "2018": 6543,
            "2019": 8765, "2020": 11896, "2021": 9876, "2022": 8654,
            "2023": 7543, "2024": 5234,
        },
        "gc_distribution": {
            "25-30": 456, "30-35": 2345, "35-40": 8765, "40-45": 12345,
            "45-50": 15678, "50-55": 14567, "55-60": 9876, "60-65": 5432,
            "65-70": 2345, "70-75": 876,
        },
        "length_distribution": {
            "0-5kb": 8765, "5-10kb": 12345, "10-25kb": 15678,
            "25-50kb": 14567, "50-100kb": 10234, "100-200kb": 6543,
            "200-500kb": 3210, ">500kb": 1018,
        },
        "mob_typing": {
            "rep_types": {
                "IncFIB": 4523, "IncFIA": 3876, "IncFII": 3654,
                "IncI1": 2345, "IncHI2": 1876, "IncN": 1654,
                "IncX4": 1432, "IncQ1": 1234, "IncR": 987,
                "IncP-1": 876, "IncL": 765, "IncC": 654,
                "IncFIC": 543, "IncB/O/K/Z": 432, "ColRNAI": 2187,
                "Col156": 1543, "Col(MG828)": 987, "rep_cluster_1": 876,
            },
            "relaxase_types": {
                "MOBF": 8765, "MOBP": 6543, "MOBQ": 3456,
                "MOBH": 1876, "MOBC": 987, "MOBV": 543,
            },
            "mpf_types": {
                "MPF_F": 7654, "MPF_T": 5432, "MPF_I": 2345,
                "MPF_G": 876, "MPF_FA": 654, "MPF_B": 432,
            },
            "mobility": {
                "conjugative": 18765, "mobilizable": 24567,
                "non-mobilizable": 29028,
            },
            "sample_size": 72360,
        },
        "correlations": {
            "sample_size": 125,
            "amr_vs_inc": {
                "BETA-LACTAM": {"IncFIB": 18, "IncFIA": 15, "IncFII": 12, "IncI1": 8, "ColRNAI": 6, "IncN": 5, "IncHI2": 4, "IncC": 3},
                "AMINOGLYCOSIDE": {"IncFIB": 14, "IncFIA": 11, "IncFII": 9, "IncI1": 7, "ColRNAI": 5, "IncN": 4, "IncHI2": 6},
                "SULFONAMIDE": {"IncFIB": 12, "IncFIA": 9, "IncFII": 8, "IncI1": 6, "IncHI2": 5, "IncN": 4},
                "TETRACYCLINE": {"IncFIB": 10, "IncFIA": 8, "IncN": 6, "IncI1": 5, "IncHI2": 4},
                "QUINOLONE": {"IncFIB": 7, "IncFIA": 6, "IncI1": 4, "IncN": 3, "IncHI2": 3},
                "TRIMETHOPRIM": {"IncFIB": 8, "IncFIA": 6, "IncI1": 5, "IncHI2": 4, "IncN": 3},
                "PHENICOL": {"IncFIB": 6, "IncHI2": 5, "IncI1": 4, "IncFIA": 3},
                "COLISTIN": {"IncI1": 4, "IncHI2": 3, "IncX4": 6, "IncFIB": 2},
                "MERCURY": {"IncFIB": 5, "IncFIA": 4, "IncHI2": 3, "IncN": 3, "IncI1": 2},
                "COPPER": {"IncFIB": 4, "IncHI2": 3, "IncFIA": 2, "IncI1": 2},
                "ARSENIC": {"IncFIB": 3, "IncHI2": 4, "IncN": 2, "IncI1": 2},
                "TELLURIUM": {"IncHI2": 5, "IncFIB": 3, "IncFIA": 2},
            },
            "amr_vs_mobility": {
                "BETA-LACTAM": {"conjugative": 28, "mobilizable": 18, "non-mobilizable": 8},
                "AMINOGLYCOSIDE": {"conjugative": 24, "mobilizable": 15, "non-mobilizable": 6},
                "SULFONAMIDE": {"conjugative": 20, "mobilizable": 12, "non-mobilizable": 5},
                "TETRACYCLINE": {"conjugative": 16, "mobilizable": 10, "non-mobilizable": 4},
                "QUINOLONE": {"conjugative": 12, "mobilizable": 8, "non-mobilizable": 3},
                "TRIMETHOPRIM": {"conjugative": 14, "mobilizable": 9, "non-mobilizable": 3},
                "PHENICOL": {"conjugative": 10, "mobilizable": 6, "non-mobilizable": 2},
                "COLISTIN": {"conjugative": 6, "mobilizable": 5, "non-mobilizable": 4},
                "MERCURY": {"conjugative": 10, "mobilizable": 4, "non-mobilizable": 1},
                "COPPER": {"conjugative": 7, "mobilizable": 3, "non-mobilizable": 1},
                "ARSENIC": {"conjugative": 6, "mobilizable": 3, "non-mobilizable": 2},
                "TELLURIUM": {"conjugative": 7, "mobilizable": 2, "non-mobilizable": 1},
            },
            "heavy_metal_genes": {
                "merA": 12, "merR": 10, "merD": 8, "merE": 7, "merP": 6,
                "arsB": 9, "arsC": 8, "arsR": 7, "arsA": 5,
                "pcoA": 6, "pcoB": 5, "pcoD": 4, "pcoR": 4, "pcoS": 3,
                "silA": 5, "silB": 4, "silC": 3, "silR": 3, "silS": 3, "silP": 2,
                "terA": 7, "terB": 6, "terC": 5, "terD": 5, "terE": 4, "terW": 3, "terZ": 3,
            },
            "heavy_metal_vs_inc": {
                "MERCURY": {"IncFIB": 5, "IncFIA": 4, "IncHI2": 3, "IncN": 3},
                "ARSENIC": {"IncHI2": 4, "IncFIB": 3, "IncN": 2, "IncI1": 2},
                "COPPER": {"IncFIB": 4, "IncHI2": 3, "IncFIA": 2},
                "COPPER/SILVER": {"IncFIB": 3, "IncHI2": 2, "IncI1": 2},
                "TELLURIUM": {"IncHI2": 5, "IncFIB": 3, "IncFIA": 2},
                "SILVER": {"IncFIB": 3, "IncHI2": 2, "IncI1": 1},
            },
            "heavy_metal_vs_mobility": {
                "MERCURY": {"conjugative": 10, "mobilizable": 4, "non-mobilizable": 1},
                "ARSENIC": {"conjugative": 6, "mobilizable": 3, "non-mobilizable": 2},
                "COPPER": {"conjugative": 7, "mobilizable": 3, "non-mobilizable": 1},
                "COPPER/SILVER": {"conjugative": 5, "mobilizable": 2, "non-mobilizable": 0},
                "TELLURIUM": {"conjugative": 7, "mobilizable": 2, "non-mobilizable": 1},
                "SILVER": {"conjugative": 4, "mobilizable": 2, "non-mobilizable": 0},
            },
            "phage_distribution": {"0": 35, "1-2": 28, "3-5": 22, "6-10": 18, ">10": 12},
            "transposase_distribution": {"0": 20, "1-3": 30, "4-8": 28, "9-15": 25, ">15": 22},
            "phage_vs_mobility": {
                "conjugative": {"mean_phage_elements": 6.8, "count": 48, "total_elements": 326},
                "mobilizable": {"mean_phage_elements": 3.2, "count": 42, "total_elements": 134},
                "non-mobilizable": {"mean_phage_elements": 1.4, "count": 35, "total_elements": 49},
            },
            "pmlst_schemes": {
                "IncF__RST": 32, "IncHI2__ST": 8, "IncI1__MLST": 7,
                "IncN__MLST": 5, "IncHI1__MLST": 3, "IncA/C__ST": 2,
            },
            "pmlst_alleles": {
                "FIA(1)": 12, "FIB(1)": 18, "FII(2)": 10, "FII(18)": 8,
                "FIA(27)": 6, "FIC(4)": 5, "FIB(20)": 4, "FII(29)": 4,
                "FII(1)": 3, "FIA(11)": 3, "FIB(10)": 3, "FIC(1)": 2,
            },
            "pmlst_vs_mobility": {
                "IncF__RST": {"conjugative": 24, "mobilizable": 6, "non-mobilizable": 2},
                "IncHI2__ST": {"conjugative": 7, "mobilizable": 1, "non-mobilizable": 0},
                "IncI1__MLST": {"conjugative": 5, "mobilizable": 2, "non-mobilizable": 0},
                "IncN__MLST": {"conjugative": 4, "mobilizable": 1, "non-mobilizable": 0},
            },
        },
    }
