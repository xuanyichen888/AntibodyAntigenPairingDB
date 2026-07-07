#!/usr/bin/env python3
"""Validate that the local project outputs are present and readable."""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

CHECKS = [
    ("master_v2", PROJECT_ROOT / "antibody_antigen_master_v2.csv"),
    ("pubmed_refs", PROJECT_ROOT / "pubmed_references.csv"),
    ("v1", PROJECT_ROOT / "outputs" / "candidates_v1.csv"),
    ("v2", PROJECT_ROOT / "outputs" / "candidates_v2.csv"),
    ("latest", PROJECT_ROOT / "outputs" / "antibody_antigen_candidates.csv"),
    ("manual_sample", PROJECT_ROOT / "outputs" / "manual_validation_sample.csv"),
    ("binder_type_review_queue", PROJECT_ROOT / "outputs" / "binder_type_review_queue.csv"),
    ("patent_search_audit", PROJECT_ROOT / "outputs" / "patent_search_audit.csv"),
    ("patent_validation_queue", PROJECT_ROOT / "outputs" / "patent_validation_queue.csv"),
    ("patent_validation_patent_summary", PROJECT_ROOT / "outputs" / "patent_validation_patent_summary.csv"),
]

REQUIRED_CANDIDATE_COLUMNS = {
    "record_id",
    "binder_name",
    "binder_sequence",
    "antigen_name",
    "target_sequence",
    "antigen_uniprot_id",
    "binder_species",
    "antigen_species",
    "normalized_KD_nM",
    "binder_type_status",
    "complex_group_id",
    "pdb_pairing_status",
    "qa_flags",
    "source_type",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> int:
    ok = True
    for label, path in CHECKS:
        if not path.exists():
            print(f"FAIL {label}: missing {path}")
            ok = False
            continue

        rows = read_rows(path)
        print(f"OK   {label}: {len(rows)} rows")
        if rows and "source_type" in rows[0]:
            print(f"     sources: {dict(Counter(r['source_type'] for r in rows))}")
        if rows and "antigen_name" in rows[0]:
            antigens = {r["antigen_name"] for r in rows if r.get("antigen_name")}
            print(f"     antigens: {len(antigens)}")

        if label in {"master_v2", "v1", "v2", "latest"}:
            fieldnames = set(rows[0].keys()) if rows else set()
            missing = sorted(REQUIRED_CANDIDATE_COLUMNS - fieldnames)
            if missing:
                print(f"FAIL {label}: missing columns {', '.join(missing)}")
                ok = False

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
