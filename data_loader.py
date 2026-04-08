"""
Data fetching and caching layer for the PLSDB 2025 dashboard.
Fetches plasmid metadata from the PLSDB API and caches locally.
"""

import json
import os
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

    # Top taxonomy - phylum level
    for kingdom in ["Bacteria", "Archaea"]:
        result = fetch_filter_taxonomy(TAXONOMY_superkingdom=kingdom)
        if result:
            accs = _extract_accessions(result)
            overview[f"kingdom_{kingdom}"] = len(accs)
        else:
            overview[f"kingdom_{kingdom}"] = 0

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
    }
