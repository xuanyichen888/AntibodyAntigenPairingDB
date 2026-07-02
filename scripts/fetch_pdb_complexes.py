#!/usr/bin/env python3
"""
Step 1a: Fetch antibody-antigen complex structures from RCSB PDB.

Searches PDB for nanobody, Fab, scFv, and antibody complexes, then extracts
polymer entity sequences and descriptions via the RCSB REST API.

Output: data/pdb/pdb_ab_ag_pairs.csv
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "data" / "pdb"
OUTPUT_PATH = OUTPUT_DIR / "pdb_ab_ag_pairs.csv"

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity"

SEARCHES = [
    "nanobody antigen complex",
    "VHH single domain antibody complex",
    "Fab antibody antigen complex immune",
    "scFv single chain variable fragment antigen",
]

AB_KEYWORDS = [
    "nanobody", "vhh", "antibody", "fab", "scfv", "heavy chain variable",
    "light chain variable", "immunoglobulin", "single domain", "variable domain",
    "sdab", "camelid", "sybody", "anti-", "anti ",
]
AG_SKIP = [
    "nanobody", "antibody", "fab fragment", "immunoglobulin",
    "light chain", "heavy chain", "constant region", "fc region", "fc fragment",
    "kappa chain", "lambda chain", "variable domain", "variable heavy",
    "variable light", "binding ig", "igm rf", "igg fc",
]

FIELDNAMES = [
    "pdb_id", "binder_name", "binder_type", "binder_sequence",
    "antigen_name", "antigen_sequence", "title", "source",
]


def rcsb_search(query_text: str, max_results: int = 500) -> list[str]:
    body = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": query_text}},
        "return_type": "entry",
        "request_options": {"return_all_hits": False, "paginate": {"start": 0, "rows": max_results}},
    }
    req = urllib.request.Request(
        SEARCH_URL, json.dumps(body).encode(), headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        result = json.loads(resp.read())
    return [r["identifier"] for r in result.get("result_set", [])]


def classify_type(desc: str) -> str:
    d = desc.lower()
    if any(kw in d for kw in ("nanobody", "vhh", "single domain", "single-domain", "sdab", "camelid", "sybody")):
        return "nanobody"
    if any(kw in d for kw in ("scfv", "single chain", "single-chain")):
        return "scFv"
    if any(kw in d for kw in ("heavy chain", " vh")):
        return "heavy chain"
    if any(kw in d for kw in ("light chain", " vl")):
        return "light chain"
    if "fab" in d:
        return "Fab"
    return "antibody"


def is_antibody(desc: str) -> bool:
    return any(kw in desc.lower() for kw in AB_KEYWORDS)


def is_antigen(desc: str) -> bool:
    d = desc.lower()
    return not any(kw in d for kw in AG_SKIP) and len(desc) > 3


def fetch_entry(pdb_id: str) -> tuple[str, list[dict]]:
    req = urllib.request.Request(f"{ENTRY_URL}/{pdb_id}")
    with urllib.request.urlopen(req, timeout=10) as resp:
        entry = json.loads(resp.read())
    title = entry.get("struct", {}).get("title", "")

    chains: list[dict] = []
    for eid in range(1, 9):
        try:
            req2 = urllib.request.Request(f"{ENTITY_URL}/{pdb_id}/{eid}")
            with urllib.request.urlopen(req2, timeout=8) as resp2:
                ent = json.loads(resp2.read())
            desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
            seq = ent.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
            ptype = ent.get("entity_poly", {}).get("type", "")
            if ptype == "polypeptide(L)" and seq and len(seq) >= 30:
                chains.append({"desc": desc, "seq": seq})
        except urllib.error.HTTPError as e:
            if e.code == 404:
                break
    return title, chains


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_ids: set[str] = set()
    for query in SEARCHES:
        try:
            ids = rcsb_search(query)
            all_ids.update(ids)
            print(f"  {query!r}: {len(ids)} structures")
        except Exception as e:
            print(f"  {query!r}: FAILED - {e}")
        time.sleep(0.3)
    print(f"Total unique PDB IDs: {len(all_ids)}")

    rows: list[dict] = []
    errors = 0
    for i, pdb_id in enumerate(sorted(all_ids)):
        try:
            title, chains = fetch_entry(pdb_id)
            ab_chains = [c for c in chains if is_antibody(c["desc"])]
            ag_chains = [c for c in chains if is_antigen(c["desc"])]
            for ab in ab_chains:
                for ag in ag_chains:
                    rows.append({
                        "pdb_id": pdb_id,
                        "binder_name": ab["desc"],
                        "binder_type": classify_type(ab["desc"]),
                        "binder_sequence": ab["seq"],
                        "antigen_name": ag["desc"],
                        "antigen_sequence": ag["seq"],
                        "title": title,
                        "source": f"PDB:{pdb_id}",
                    })
        except Exception:
            errors += 1
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(all_ids)}, pairs: {len(rows)}, errors: {errors}")
        time.sleep(0.15)

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {OUTPUT_PATH} ({len(rows)} rows, {errors} errors)")
    types = Counter(r["binder_type"] for r in rows)
    print(f"Binder types: {dict(types)}")
    print(f"Unique antigens: {len(set(r['antigen_name'] for r in rows))}")


if __name__ == "__main__":
    main()
