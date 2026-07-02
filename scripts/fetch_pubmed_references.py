#!/usr/bin/env python3
"""
Fetch PubMed references for antibody-antigen candidate targets.

PubMed is used here as a literature evidence layer, not as a direct
sequence-pair source. Abstract metadata rarely contains paired binder and
antigen amino-acid sequences, so this script writes reference rows that can be
used during manual curation and in the Streamlit app.

Default input:  outputs/antibody_antigen_candidates.csv
Default output: pubmed_references.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import ssl
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "outputs" / "antibody_antigen_candidates.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "pubmed_references.csv"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

try:
    import certifi
except ImportError:  # pragma: no cover
    SSL_CONTEXT = ssl.create_default_context()
else:
    SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())

FIELDNAMES = [
    "antigen_name",
    "query",
    "pmid",
    "title",
    "journal",
    "pubdate",
    "doi",
    "authors",
    "source_reference",
    "evidence_text",
    "url",
]

ANTIBODY_TERMS = [
    "antibody",
    "antibodies",
    "nanobody",
    "single-domain antibody",
    "scFv",
    "Fab",
    "neutralizing antibody",
]


def clean(value: str | None) -> str:
    return (value or "").strip()


def entrez_get(endpoint: str, params: dict[str, Any], timeout: int = 30) -> dict[str, Any]:
    merged = dict(params)
    api_key = os.environ.get("NCBI_API_KEY") or os.environ.get("NCBI_EUTILS_API_KEY")
    email = os.environ.get("NCBI_EMAIL")
    tool = os.environ.get("NCBI_TOOL", "antigen_antibody_pairing_pubmed")
    if api_key:
        merged["api_key"] = api_key
    if email:
        merged["email"] = email
    if tool:
        merged["tool"] = tool

    query = urllib.parse.urlencode(merged)
    url = f"{EUTILS_BASE}/{endpoint}.fcgi?{query}"
    req = urllib.request.Request(url, headers={"User-Agent": tool})
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CONTEXT) as resp:
        return json.loads(resp.read())


def load_targets(path: Path, max_targets: int, min_count: int) -> list[str]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    counts = Counter(clean(row.get("antigen_name")) for row in rows)
    counts.pop("", None)
    targets = [name for name, count in counts.most_common() if count >= min_count]
    if max_targets > 0:
        targets = targets[:max_targets]
    return targets


def build_query(antigen_name: str) -> str:
    antigen = antigen_name.replace('"', "")
    antibody_clause = " OR ".join(f'"{term}"[Title/Abstract]' for term in ANTIBODY_TERMS)
    antigen_clause = f'"{antigen}"[Title/Abstract] OR "{antigen}"[MeSH Terms]'
    return f"({antigen_clause}) AND ({antibody_clause})"


def search_pubmed(query: str, retmax: int) -> list[str]:
    data = entrez_get(
        "esearch",
        {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": retmax,
            "sort": "relevance",
        },
    )
    return data.get("esearchresult", {}).get("idlist", [])


def summarize_pmids(pmids: list[str]) -> list[dict[str, Any]]:
    if not pmids:
        return []
    data = entrez_get(
        "esummary",
        {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
        },
    )
    result = data.get("result", {})
    return [result[pmid] for pmid in result.get("uids", []) if pmid in result]


def article_id(record: dict[str, Any], idtype: str) -> str:
    for item in record.get("articleids", []):
        if item.get("idtype") == idtype:
            return clean(item.get("value"))
    return ""


def authors(record: dict[str, Any], limit: int = 6) -> str:
    names = [clean(author.get("name")) for author in record.get("authors", [])]
    names = [name for name in names if name]
    if len(names) > limit:
        return "; ".join(names[:limit]) + "; et al."
    return "; ".join(names)


def record_to_row(antigen_name: str, query: str, record: dict[str, Any]) -> dict[str, str]:
    pmid = clean(record.get("uid"))
    doi = article_id(record, "doi")
    title = clean(record.get("title"))
    journal = clean(record.get("fulljournalname")) or clean(record.get("source"))
    pubdate = clean(record.get("pubdate"))
    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
    evidence_bits = [title, journal, pubdate]
    return {
        "antigen_name": antigen_name,
        "query": query,
        "pmid": pmid,
        "title": title,
        "journal": journal,
        "pubdate": pubdate,
        "doi": doi,
        "authors": authors(record),
        "source_reference": f"PMID:{pmid}" if pmid else "",
        "evidence_text": " | ".join(bit for bit in evidence_bits if bit),
        "url": url,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--targets", default="", help="Comma-separated antigen names. Overrides --input.")
    parser.add_argument("--max-targets", type=int, default=30, help="Top targets to query. Use 0 for all.")
    parser.add_argument("--min-count", type=int, default=2, help="Minimum candidate rows per target.")
    parser.add_argument("--retmax", type=int, default=5, help="PubMed records per target.")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between Entrez calls.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.targets.strip():
        targets = [clean(t) for t in args.targets.split(",") if clean(t)]
    else:
        targets = load_targets(args.input, args.max_targets, args.min_count)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    failures = 0

    print(f"Querying PubMed for {len(targets)} targets...")
    for target in targets:
        query = build_query(target)
        try:
            pmids = search_pubmed(query, args.retmax)
            time.sleep(args.sleep)
            records = summarize_pmids(pmids)
        except Exception as exc:  # noqa: BLE001 - keep batch runs moving.
            failures += 1
            print(f"  {target}: failed - {exc}")
            time.sleep(args.sleep)
            continue

        for record in records:
            row = record_to_row(target, query, record)
            key = (row["antigen_name"], row["pmid"])
            if row["pmid"] and key not in seen:
                rows.append(row)
                seen.add(key)
        print(f"  {target}: {len(records)} references")
        time.sleep(args.sleep)

    if targets and failures == len(targets) and not rows:
        raise SystemExit("All PubMed requests failed; output was not overwritten.")

    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {args.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
