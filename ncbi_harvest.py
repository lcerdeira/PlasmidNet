#!/usr/bin/env python3
"""
Harvest complete plasmid metadata from NCBI Nucleotide.
Downloads accession, length, topology, organism, dates, and BioSample IDs
for all complete plasmids not already in PLSDB.

Usage:
    python ncbi_harvest.py          # full harvest (~30 min)
    python ncbi_harvest.py --test   # test with 100 records

Outputs data/ncbi_plasmids.csv which is then imported by build_db.py.
"""

import csv
import os
import sqlite3
import sys
import time
import xml.etree.ElementTree as ET

import requests

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

DB_PATH = os.path.join(os.path.dirname(__file__), "data", "plsdb.db")
OUT_PATH = os.path.join(os.path.dirname(__file__), "data", "raw", "ncbi_plasmids.csv")
BATCH_SIZE = 500
API_KEY = os.environ.get("NCBI_API_KEY", "")  # optional, increases rate limit


def esearch_all_plasmids():
    """Search NCBI for all complete plasmid sequences. Returns (count, webenv, query_key)."""
    params = {
        "db": "nuccore",
        "term": (
            "plasmid[title] AND complete[title] "
            "AND biomol_genomic[prop] "
            'AND ("1000"[SLEN] : "10000000"[SLEN])'  # 1kb - 10Mb
        ),
        "retmax": 0,
        "usehistory": "y",
    }
    if API_KEY:
        params["api_key"] = API_KEY

    resp = requests.get(NCBI_ESEARCH, params=params, timeout=60)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    count = int(root.findtext("Count", "0"))
    webenv = root.findtext("WebEnv", "")
    query_key = root.findtext("QueryKey", "1")
    return count, webenv, query_key


def efetch_batch(webenv, query_key, retstart, retmax=BATCH_SIZE):
    """Fetch a batch of document summaries from NCBI history server."""
    params = {
        "db": "nuccore",
        "query_key": query_key,
        "WebEnv": webenv,
        "retstart": retstart,
        "retmax": retmax,
        "rettype": "docsum",
        "retmode": "xml",
    }
    if API_KEY:
        params["api_key"] = API_KEY

    for attempt in range(3):
        try:
            resp = requests.get(NCBI_ESUMMARY, params=params, timeout=120)
            resp.raise_for_status()
            return resp.text
        except Exception as e:
            if attempt < 2:
                print(f"    Retry {attempt + 1} after error: {e}")
                time.sleep(5 * (attempt + 1))
            else:
                raise


def parse_docsum(xml_text):
    """Parse NCBI DocSum XML into a list of plasmid records."""
    records = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return records

    for doc in root.findall(".//DocSum"):
        uid = doc.findtext("Id", "")
        items = {}
        for item in doc.findall("Item"):
            name = item.get("Name", "")
            val = item.text or ""
            items[name] = val
            # Handle sub-items (e.g., BioSample in extra)
            for sub in item.findall("Item"):
                items[sub.get("Name", "")] = sub.text or ""

        acc = items.get("AccessionVersion", items.get("Caption", ""))
        title = items.get("Title", "")
        length = int(items.get("Length", 0) or 0)
        create_date = items.get("CreateDate", "")
        update_date = items.get("UpdateDate", "")
        taxid = items.get("TaxId", "")
        organism = items.get("Organism", "")
        # Extract organism from title if not in Organism field
        if not organism and " strain " in title:
            organism = title.split(" strain ")[0].strip()
        elif not organism and " plasmid " in title:
            organism = title.split(" plasmid ")[0].strip()
        # Topology from extra field or title
        extra = items.get("Extra", "")
        if "topology=circular" in extra.lower():
            topology = "circular"
        elif "topology=linear" in extra.lower():
            topology = "linear"
        else:
            topology = "circular" if "circular" in title.lower() else "linear"

        # Skip non-plasmid entries
        if length < 1000 or length > 10_000_000:
            continue

        records.append({
            "NCBI_UID": uid,
            "accession": acc,
            "description": title,
            "length": length,
            "topology": topology,
            "create_date": create_date,
            "update_date": update_date,
            "taxid": taxid,
            "organism": organism,
        })

    return records


def get_existing_accessions():
    """Get set of accessions already in PLSDB."""
    if not os.path.exists(DB_PATH):
        return set()
    conn = sqlite3.connect(DB_PATH)
    rows = conn.execute("SELECT NUCCORE_ACC FROM nuccore").fetchall()
    conn.close()
    # Store base accession (without version) for flexible matching
    existing = set()
    for (acc,) in rows:
        existing.add(acc)
        existing.add(acc.split(".")[0])  # also match without version
    return existing


def main():
    test_mode = "--test" in sys.argv
    max_records = 200 if test_mode else None

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

    print("=== NCBI Complete Plasmid Harvester ===")

    # Step 1: Search
    print("Searching NCBI for complete plasmids...")
    count, webenv, query_key = esearch_all_plasmids()
    print(f"  Found {count:,} complete plasmids in NCBI")

    if max_records:
        count = min(count, max_records)
        print(f"  Test mode: limiting to {count}")

    # Step 2: Get existing PLSDB accessions
    existing = get_existing_accessions()
    print(f"  PLSDB has {len(existing):,} existing accessions")

    # Step 3: Fetch in batches
    print(f"  Fetching metadata in batches of {BATCH_SIZE}...")
    all_records = []
    new_count = 0
    skipped = 0

    for start in range(0, count, BATCH_SIZE):
        batch_num = start // BATCH_SIZE + 1
        total_batches = (count + BATCH_SIZE - 1) // BATCH_SIZE
        print(f"  Batch {batch_num}/{total_batches} "
              f"(records {start + 1}-{min(start + BATCH_SIZE, count)})...")

        xml_text = efetch_batch(webenv, query_key, start, BATCH_SIZE)
        records = parse_docsum(xml_text)

        for rec in records:
            acc_base = rec["accession"].split(".")[0]
            if rec["accession"] in existing or acc_base in existing:
                skipped += 1
                continue
            all_records.append(rec)
            new_count += 1

        # Rate limiting: 3 requests/sec without key, 10/sec with key
        time.sleep(0.35 if not API_KEY else 0.1)

    print(f"\n  Total fetched: {new_count + skipped:,}")
    print(f"  Already in PLSDB: {skipped:,}")
    print(f"  New plasmids: {new_count:,}")

    # Step 4: Write CSV
    if all_records:
        fieldnames = ["NCBI_UID", "accession", "description", "length",
                      "topology", "create_date", "update_date",
                      "taxid", "organism"]
        with open(OUT_PATH, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_records)
        print(f"\n  Wrote {len(all_records):,} records to {OUT_PATH}")
    else:
        print("\n  No new plasmids found")

    # Step 5: Summary stats
    if all_records:
        organisms = {}
        years = {}
        for r in all_records:
            org = r["organism"].split()[0] if r["organism"] else "Unknown"
            organisms[org] = organisms.get(org, 0) + 1
            year = r["create_date"][:4] if r["create_date"] else "?"
            years[year] = years.get(year, 0) + 1

        print("\n  Top 10 organisms:")
        for org, cnt in sorted(organisms.items(), key=lambda x: -x[1])[:10]:
            print(f"    {org:30s} {cnt:6,}")

        print("\n  By year:")
        for year, cnt in sorted(years.items()):
            if year.isdigit() and int(year) >= 2020:
                print(f"    {year}: {cnt:,}")

    print("\nDone!")


if __name__ == "__main__":
    main()
