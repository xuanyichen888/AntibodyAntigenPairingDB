#!/usr/bin/env python3
"""
Step 1b: Fetch patent-deposited antibody sequences from NCBI GenBank PAT division.

Searches NCBI Entrez protein database for antibody variable regions targeting
specific antigens, then fetches canonical antigen sequences from UniProt.

Output: data/patent/patent_ab_ag_pairs.csv
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "data" / "patent"
OUTPUT_PATH = OUTPUT_DIR / "patent_ab_ag_pairs.csv"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"

TARGETS = {
    "PD-1": {"query": "anti-PD-1 VH OR anti-PD-1 VHH OR anti-PD-1 scFv", "uniprot": "Q15116"},
    "PD-L1": {"query": "anti-PD-L1 VH OR anti-PD-L1 VHH OR pembrolizumab", "uniprot": "Q9NZQ7"},
    "EGFR": {"query": "anti-EGFR VH OR cetuximab heavy chain", "uniprot": "P00533"},
    "TNF-alpha": {"query": "anti-TNF VH OR adalimumab heavy chain", "uniprot": "P01375"},
    "HER2": {"query": "anti-HER2 VH OR trastuzumab heavy chain", "uniprot": "P04626"},
    "VEGF": {"query": "anti-VEGF VH OR bevacizumab heavy chain", "uniprot": "P15692"},
    "CD20": {"query": "anti-CD20 VH OR rituximab heavy chain", "uniprot": "P11836"},
    "CD19": {"query": "anti-CD19 VH OR anti-CD19 scFv", "uniprot": "P15391"},
    "CD3": {"query": "anti-CD3 VH OR anti-CD3 scFv", "uniprot": "P07766"},
    "RSV": {"query": "anti-RSV VH OR palivizumab heavy chain", "uniprot": "P03420"},
    "SARS-CoV-2 Spike": {"query": "anti-SARS-CoV-2 VH OR anti-spike antibody", "uniprot": "P0DTC2"},
    "TIGIT": {"query": "anti-TIGIT VH OR anti-TIGIT VHH", "uniprot": "Q495A1"},
    "LAG-3": {"query": "anti-LAG-3 VH OR anti-LAG3 antibody", "uniprot": "P18627"},
    "TIM-3": {"query": "anti-TIM-3 VH OR anti-TIM3 antibody", "uniprot": "Q8TDQ0"},
    "CD47": {"query": "anti-CD47 VH OR anti-CD47 VHH", "uniprot": "Q08722"},
    "CTLA-4": {"query": "anti-CTLA-4 VH OR ipilimumab heavy chain", "uniprot": "P16410"},
    "NGF": {"query": "anti-NGF VH OR anti-NGF antibody", "uniprot": "P01138"},
}

FIELDNAMES = [
    "ncbi_accession", "binder_name", "binder_type", "binder_sequence",
    "antigen_name", "antigen_sequence", "patent_id", "source",
]

AA_CHARS = set("ACDEFGHIKLMNPQRSTVWYXBZU")


def parse_accession(header: str) -> str:
    token = header.split()[0]
    parts = token.split("|")
    for part in parts:
        if re.fullmatch(r"[A-Z]{2,4}_?\d+(?:\.\d+)?", part):
            return part
    for part in parts:
        if part and part.lower() not in {"gb", "emb", "dbj", "pat", "ref", "sp", "tr", "us"}:
            return part
    return token


def extract_patent_id(title: str) -> str:
    match = re.search(r"\b(?:US|WO|EP)\s*\d+[A-Z0-9 ]{0,8}", title, re.IGNORECASE)
    if match:
        return re.sub(r"\s+", " ", match.group(0)).strip(" .;,")
    match = re.search(r"Patent:\s*([^,\]\s]+)", title, re.IGNORECASE)
    return match.group(1) if match else ""


def classify_binder(title: str) -> str:
    t = title.lower()
    if any(kw in t for kw in ("nanobody", "vhh", "single domain", "single-domain", "sdab")):
        return "nanobody"
    if any(kw in t for kw in ("scfv", "single chain", "single-chain")):
        return "scFv"
    if any(kw in t for kw in ("heavy chain", " vh ")):
        return "heavy chain"
    if any(kw in t for kw in ("light chain", " vl ")):
        return "light chain"
    if "fab" in t:
        return "Fab"
    return "antibody"


def ncbi_search(query: str, retmax: int = 50) -> list[str]:
    params = urllib.parse.urlencode({
        "db": "protein",
        "term": f'({query}) AND "PAT"[Division]',
        "retmax": retmax,
        "retmode": "json",
    })
    url = f"{ESEARCH_URL}?{params}"
    with urllib.request.urlopen(url, timeout=15) as resp:
        data = json.loads(resp.read())
    return data.get("esearchresult", {}).get("idlist", [])


def ncbi_fetch_fasta(gis: list[str]) -> list[dict]:
    if not gis:
        return []
    params = urllib.parse.urlencode({
        "db": "protein",
        "id": ",".join(gis),
        "rettype": "fasta",
        "retmode": "text",
    })
    url = f"{EFETCH_URL}?{params}"
    with urllib.request.urlopen(url, timeout=30) as resp:
        text = resp.read().decode()

    entries = []
    for block in text.strip().split(">"):
        block = block.strip()
        if not block:
            continue
        lines = block.split("\n")
        header = lines[0]
        seq = "".join(line.strip() for line in lines[1:])
        accession = parse_accession(header)
        entries.append({"accession": accession, "title": header, "sequence": seq})
    return entries


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    url = f"{UNIPROT_URL}/{uniprot_id}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            text = resp.read().decode()
        lines = text.strip().split("\n")
        return "".join(line.strip() for line in lines[1:])
    except Exception:
        return ""


def is_antibody_sequence(seq: str, title: str) -> bool:
    if len(seq) < 80 or len(seq) > 600:
        return False
    if not all(c in AA_CHARS for c in seq.upper()):
        return False
    t = title.lower()
    return any(kw in t for kw in (
        "antibod", "immunoglobulin", "variable", "vh", "vl", "vhh",
        "nanobody", "scfv", "fab", "cdr", "heavy chain", "light chain",
        "anti-",
    ))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    antigen_seqs: dict[str, str] = {}
    print("Fetching antigen sequences from UniProt...")
    for name, info in TARGETS.items():
        seq = fetch_uniprot_sequence(info["uniprot"])
        if seq:
            antigen_seqs[name] = seq
            print(f"  {name}: {len(seq)} aa")
        else:
            print(f"  {name}: FAILED")
        time.sleep(0.3)

    rows: list[dict] = []
    for target_name, info in TARGETS.items():
        ag_seq = antigen_seqs.get(target_name, "")
        if not ag_seq:
            continue
        try:
            gis = ncbi_search(info["query"])
            print(f"  {target_name}: {len(gis)} hits")
        except Exception as e:
            print(f"  {target_name}: search failed - {e}")
            continue
        time.sleep(1.0)

        if not gis:
            continue
        try:
            entries = ncbi_fetch_fasta(gis)
        except Exception as e:
            print(f"  {target_name}: fetch failed - {e}")
            continue
        time.sleep(1.0)

        for entry in entries:
            if is_antibody_sequence(entry["sequence"], entry["title"]):
                patent_id = extract_patent_id(entry["title"])
                rows.append({
                    "ncbi_accession": entry["accession"],
                    "binder_name": entry["title"].split("|")[-1].strip()[:200] if "|" in entry["title"] else entry["title"][:200],
                    "binder_type": classify_binder(entry["title"]),
                    "binder_sequence": entry["sequence"],
                    "antigen_name": target_name,
                    "antigen_sequence": ag_seq,
                    "patent_id": patent_id,
                    "source": f"GenBank:{entry['accession']}",
                })

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {OUTPUT_PATH} ({len(rows)} rows)")
    types = Counter(r["binder_type"] for r in rows)
    targets = Counter(r["antigen_name"] for r in rows)
    print(f"Binder types: {dict(types)}")
    print(f"Targets: {dict(targets)}")


if __name__ == "__main__":
    main()
