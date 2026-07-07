#!/usr/bin/env python3
"""Build patent search audit and manual validation queue artifacts.

These files bridge the gap between the automated NCBI Protein PAT sequence
retrieval and the manual USPTO / Google Patents validation step. They do not
claim that USPTO or Google Patents have been parsed automatically; they create
review-ready links, priority labels, and blank curation fields.
"""

from __future__ import annotations

import csv
import re
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

from fetch_patent_sequences import AB_FILTER, TARGETS


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MASTER_PATH = PROJECT_ROOT / "antibody_antigen_master_v2.csv"
PATENT_PATH = PROJECT_ROOT / "data" / "patent" / "patent_ab_ag_pairs.csv"
SEARCH_STRATEGY_PATH = PROJECT_ROOT / "docs" / "search_strategy.md"
OUTPUT_DIR = PROJECT_ROOT / "outputs"

SEARCH_AUDIT_PATH = OUTPUT_DIR / "patent_search_audit.csv"
VALIDATION_QUEUE_PATH = OUTPUT_DIR / "patent_validation_queue.csv"
PATENT_SUMMARY_PATH = OUTPUT_DIR / "patent_validation_patent_summary.csv"

RETMAX = 200
GOOGLE_AB_QUERY = "antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody OR Fab OR CDR"
PATENT_NUMBER_RE = re.compile(
    r"\b(?:US|WO|EP)\s*\d{4,12}\s*(?:[A-Z]\d?|[A-Z]{1,2}\d?)?\b",
    re.IGNORECASE,
)
SEQ_LABEL_RE = re.compile(r"\bSequence\s+(\d+)\b", re.IGNORECASE)

TARGET_CATEGORIES = {
    "PD-1": "Checkpoint", "PD-L1": "Checkpoint", "CTLA-4": "Checkpoint",
    "TIGIT": "Checkpoint", "LAG-3": "Checkpoint", "TIM-3": "Checkpoint",
    "VISTA": "Checkpoint", "BTLA": "Checkpoint", "CD47": "Checkpoint",
    "SIRPa": "Checkpoint",
    "HER2": "Oncology", "EGFR": "Oncology", "HER3": "Oncology",
    "CD20": "Oncology", "CD19": "Oncology", "CD3": "Oncology",
    "CD38": "Oncology", "BCMA": "Oncology", "GPRC5D": "Oncology",
    "FcRH5": "Oncology", "CD33": "Oncology", "CD22": "Oncology",
    "CD30": "Oncology", "TROP2": "Oncology", "Nectin-4": "Oncology",
    "CLDN18.2": "Oncology", "DLL3": "Oncology", "Mesothelin": "Oncology",
    "MET": "Oncology", "PSMA": "Oncology", "FAP": "Oncology",
    "GPC3": "Oncology", "B7-H3": "Oncology", "EpCAM": "Oncology",
    "TNF-alpha": "Cytokine", "IL-6": "Cytokine", "IL-6R": "Cytokine",
    "IL-17A": "Cytokine", "IL-17RA": "Cytokine", "IL-4Ra": "Cytokine",
    "IL-13": "Cytokine", "IL-5": "Cytokine", "IL-5Ra": "Cytokine",
    "IL-23": "Cytokine", "IL-12/23 p40": "Cytokine", "IL-33": "Cytokine",
    "IL-31": "Cytokine", "TSLP": "Cytokine", "IgE": "Cytokine",
    "VEGF": "Soluble", "VEGFR2": "Soluble", "NGF": "Soluble",
    "CGRP": "Soluble", "CGRPr": "Soluble", "PCSK9": "Soluble",
    "C5": "Soluble", "RANKL": "Soluble", "Ang2": "Soluble", "BAFF": "Soluble",
    "RSV F": "Infectious", "SARS-CoV-2 Spike": "Infectious",
    "Influenza HA": "Infectious", "HIV gp120": "Infectious", "Ebola GP": "Infectious",
}

SEARCH_AUDIT_COLUMNS = [
    "target_name", "category", "target_query", "uniprot_id", "full_ncbi_query",
    "retmax", "ncbi_hits_reported", "hit_count_capped_at_retmax",
    "patent_intermediate_rows", "final_master_rows", "unique_patents_final",
    "unique_accessions_final", "manual_ncbi_url", "google_patents_query",
    "google_patents_url", "uspto_manual_query", "manual_validation_status",
]

VALIDATION_QUEUE_COLUMNS = [
    "priority", "record_id", "antigen_name", "ncbi_accession",
    "patent_number", "patent_number_normalized", "seq_id_no_guess",
    "binder_name", "binder_type", "binder_type_status", "binder_sequence_length",
    "source_reference", "source_dataset", "qa_flags", "manual_ncbi_url",
    "google_patents_url", "uspto_search_url", "validation_questions",
    "manual_patent_family_verified", "manual_seq_id_no_verified",
    "manual_chain_role", "manual_claim_context", "manual_pairing_notes",
    "manual_curator", "manual_date",
]

PATENT_SUMMARY_COLUMNS = [
    "priority", "patent_number", "patent_number_normalized", "row_count",
    "target_count", "targets", "unresolved_binder_type_rows", "sequence_ids_seen",
    "ncbi_accessions", "google_patents_url", "uspto_search_url",
    "manual_validation_status",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_patent_number(value: str) -> str:
    match = PATENT_NUMBER_RE.search(value or "")
    if not match:
        return ""
    raw = re.sub(r"\s+", "", match.group(0).upper())
    return raw


def patent_display(value: str) -> str:
    normalized = normalize_patent_number(value)
    if not normalized:
        return ""
    match = re.match(r"^(US|WO|EP)(\d+)([A-Z]\d?|[A-Z]{1,2}\d?)?$", normalized)
    if not match:
        return normalized
    country, number, suffix = match.groups()
    return " ".join(part for part in [country, number, suffix or ""] if part)


def seq_id_guess(text: str) -> str:
    match = SEQ_LABEL_RE.search(text or "")
    return match.group(1) if match else ""


def encoded_url(base: str, query: str) -> str:
    return f"{base}{urllib.parse.quote_plus(query)}"


def ncbi_url(query: str) -> str:
    return encoded_url("https://www.ncbi.nlm.nih.gov/protein/?term=", query)


def google_patents_url(query: str) -> str:
    return encoded_url("https://patents.google.com/?q=", query)


def uspto_search_url() -> str:
    return "https://www.uspto.gov/patents/search/patent-public-search"


def parse_search_strategy_counts() -> dict[str, tuple[str, str]]:
    counts: dict[str, tuple[str, str]] = {}
    if not SEARCH_STRATEGY_PATH.exists():
        return counts
    for line in SEARCH_STRATEGY_PATH.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|") or "Target" in line or "---" in line:
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if len(cells) < 6:
            continue
        target, _category, _synonyms, _uniprot, hits, kept = cells[:6]
        if target and target != "":
            counts[target] = (hits.replace(",", ""), kept.replace(",", ""))
    return counts


def priority_for_row(row: dict[str, str]) -> str:
    status = row.get("binder_type_status", "")
    btype = row.get("binder_type", "")
    flags = row.get("qa_flags", "")
    length = len(row.get("binder_sequence", ""))
    if status == "non_antibody_candidate":
        return "P0_exclude_or_reclassify"
    if status == "needs_manual_review":
        return "P1_binder_type_and_seq_id"
    if btype in {"heavy chain", "light chain"}:
        return "P1_chain_pairing"
    if "sequence_contains_x" in flags or length < 80 or length > 600:
        return "P2_sequence_quality"
    return "P3_patent_context"


def validation_questions(row: dict[str, str], seq_guess: str) -> str:
    questions = [
        "confirm patent family",
        "confirm SEQ ID NO" if seq_guess else "find SEQ ID NO",
        "assign chain role",
        "check claimed/example/background context",
    ]
    if row.get("binder_type_status") == "needs_manual_review":
        questions.insert(2, "resolve binder type")
    return "; ".join(questions)


def build_search_audit(master_rows: list[dict[str, str]], patent_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    doc_counts = parse_search_strategy_counts()
    master_patent = [r for r in master_rows if r.get("source_type") == "NCBI Patent"]
    final_by_target = Counter(r.get("antigen_name", "") for r in master_patent)
    accessions_by_target: dict[str, set[str]] = defaultdict(set)
    patents_by_target: dict[str, set[str]] = defaultdict(set)
    intermediate_by_target = Counter(r.get("antigen_name", "") for r in patent_rows)

    for row in master_patent:
        target = row.get("antigen_name", "")
        accessions_by_target[target].add(row.get("source_id", ""))
        patent = normalize_patent_number(row.get("source_reference", ""))
        if patent:
            patents_by_target[target].add(patent)

    rows = []
    for target, info in TARGETS.items():
        target_query = info["query"]
        full_query = f'({target_query}) AND {AB_FILTER} AND "PAT"[Division]'
        google_query = f"({target_query}) ({GOOGLE_AB_QUERY})"
        hits, kept = doc_counts.get(target, ("", ""))
        capped = "yes" if hits.isdigit() and int(hits) >= RETMAX else "no"
        rows.append({
            "target_name": target,
            "category": TARGET_CATEGORIES.get(target, "Other"),
            "target_query": target_query,
            "uniprot_id": info["uniprot"],
            "full_ncbi_query": full_query,
            "retmax": str(RETMAX),
            "ncbi_hits_reported": hits,
            "hit_count_capped_at_retmax": capped,
            "patent_intermediate_rows": str(intermediate_by_target.get(target, 0)),
            "final_master_rows": str(final_by_target.get(target, 0)),
            "unique_patents_final": str(len(patents_by_target.get(target, set()))),
            "unique_accessions_final": str(len(accessions_by_target.get(target, set()))),
            "manual_ncbi_url": ncbi_url(full_query),
            "google_patents_query": google_query,
            "google_patents_url": google_patents_url(google_query),
            "uspto_manual_query": f'"{target}" antibody variable region',
            "manual_validation_status": "not_started",
        })
    return rows


def build_validation_queue(master_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    queue = []
    for row in master_rows:
        if row.get("source_type") != "NCBI Patent":
            continue
        patent_display_value = patent_display(row.get("source_reference", ""))
        patent_normalized = normalize_patent_number(row.get("source_reference", ""))
        seq_guess = seq_id_guess(row.get("binder_name", ""))
        search_query = patent_display_value or row.get("source_id", "")
        queue.append({
            "priority": priority_for_row(row),
            "record_id": row.get("record_id", ""),
            "antigen_name": row.get("antigen_name", ""),
            "ncbi_accession": row.get("source_id", ""),
            "patent_number": patent_display_value,
            "patent_number_normalized": patent_normalized,
            "seq_id_no_guess": seq_guess,
            "binder_name": row.get("binder_name", ""),
            "binder_type": row.get("binder_type", ""),
            "binder_type_status": row.get("binder_type_status", ""),
            "binder_sequence_length": str(len(row.get("binder_sequence", ""))),
            "source_reference": row.get("source_reference", ""),
            "source_dataset": row.get("source_dataset", ""),
            "qa_flags": row.get("qa_flags", ""),
            "manual_ncbi_url": f"https://www.ncbi.nlm.nih.gov/protein/{row.get('source_id', '')}",
            "google_patents_url": google_patents_url(search_query),
            "uspto_search_url": uspto_search_url(),
            "validation_questions": validation_questions(row, seq_guess),
            "manual_patent_family_verified": "",
            "manual_seq_id_no_verified": "",
            "manual_chain_role": "",
            "manual_claim_context": "",
            "manual_pairing_notes": "",
            "manual_curator": "",
            "manual_date": "",
        })
    queue.sort(key=lambda r: (r["priority"], r["antigen_name"], r["patent_number_normalized"], r["record_id"]))
    return queue


def build_patent_summary(queue_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in queue_rows:
        key = row["patent_number_normalized"] or row["ncbi_accession"]
        groups[key].append(row)

    rows = []
    priority_rank = {
        "P0_exclude_or_reclassify": 0,
        "P1_binder_type_and_seq_id": 1,
        "P1_chain_pairing": 1,
        "P2_sequence_quality": 2,
        "P3_patent_context": 3,
    }
    for _key, group in groups.items():
        priorities = sorted({r["priority"] for r in group}, key=lambda p: priority_rank.get(p, 9))
        targets = sorted({r["antigen_name"] for r in group if r["antigen_name"]})
        seq_ids = sorted({r["seq_id_no_guess"] for r in group if r["seq_id_no_guess"]}, key=lambda s: int(s))
        accessions = sorted({r["ncbi_accession"] for r in group if r["ncbi_accession"]})
        patent_number = group[0]["patent_number"]
        patent_normalized = group[0]["patent_number_normalized"]
        rows.append({
            "priority": priorities[0] if priorities else "",
            "patent_number": patent_number,
            "patent_number_normalized": patent_normalized,
            "row_count": str(len(group)),
            "target_count": str(len(targets)),
            "targets": "; ".join(targets),
            "unresolved_binder_type_rows": str(sum(1 for r in group if r["binder_type_status"] == "needs_manual_review")),
            "sequence_ids_seen": "; ".join(seq_ids[:80]),
            "ncbi_accessions": "; ".join(accessions[:80]),
            "google_patents_url": google_patents_url(patent_number or accessions[0] if accessions else ""),
            "uspto_search_url": uspto_search_url(),
            "manual_validation_status": "not_started",
        })
    rows.sort(key=lambda r: (
        priority_rank.get(r["priority"], 9),
        -int(r["row_count"]),
        r["patent_number_normalized"],
    ))
    return rows


def main() -> None:
    master_rows = read_csv(MASTER_PATH)
    patent_rows = read_csv(PATENT_PATH)

    search_audit = build_search_audit(master_rows, patent_rows)
    validation_queue = build_validation_queue(master_rows)
    patent_summary = build_patent_summary(validation_queue)

    write_csv(SEARCH_AUDIT_PATH, SEARCH_AUDIT_COLUMNS, search_audit)
    write_csv(VALIDATION_QUEUE_PATH, VALIDATION_QUEUE_COLUMNS, validation_queue)
    write_csv(PATENT_SUMMARY_PATH, PATENT_SUMMARY_COLUMNS, patent_summary)

    print(f"Wrote {SEARCH_AUDIT_PATH} ({len(search_audit)} rows)")
    print(f"Wrote {VALIDATION_QUEUE_PATH} ({len(validation_queue)} rows)")
    print(f"Wrote {PATENT_SUMMARY_PATH} ({len(patent_summary)} rows)")
    print(f"Priorities: {dict(Counter(r['priority'] for r in validation_queue))}")


if __name__ == "__main__":
    main()
