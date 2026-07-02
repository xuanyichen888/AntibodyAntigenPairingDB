#!/usr/bin/env python3
"""
Step 3: Merge all sources into a unified candidates table.

Reads:
  - data/source_pilots/antibody_antigen_pilot.csv  (ProteinBase seed)
  - data/pdb/pdb_ab_ag_pairs.csv                   (PDB complexes)
  - data/patent/patent_ab_ag_pairs.csv              (patent sequences)

Writes:
  - outputs/antibody_antigen_candidates.csv         (merged, deduplicated)
  - outputs/candidates_v1.csv                       (seed only, 300 rows)
  - outputs/candidates_v2.csv                       (all sources merged)
  - outputs/manual_validation_sample.csv            (stratified sample)
  - outputs/validation_memo.md                      (updated stats)
"""

from __future__ import annotations

import csv
import hashlib
import random
import re
from datetime import date
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SEED_PATH = PROJECT_ROOT / "data" / "source_pilots" / "antibody_antigen_pilot.csv"
PDB_PATH = PROJECT_ROOT / "data" / "pdb" / "pdb_ab_ag_pairs.csv"
PATENT_PATH = PROJECT_ROOT / "data" / "patent" / "patent_ab_ag_pairs.csv"
LITERATURE_PATH = PROJECT_ROOT / "pubmed_references.csv"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
MASTER_PATH = PROJECT_ROOT / "antibody_antigen_master_v2.csv"

try:
    from normalize_antigens import normalize_antigen
except ImportError:  # pragma: no cover - keeps this script runnable from odd cwd setups.
    def normalize_antigen(name: str) -> str:
        return name

CANDIDATE_COLUMNS = [
    "record_id", "binder_name", "binder_type", "binder_sequence",
    "antigen_name", "target_name", "target_sequence",
    "interaction_label", "affinity_value", "affinity_unit",
    "evidence_text", "source_reference", "confidence", "needs_review",
    "source_type", "source_dataset", "source_id", "source_confidence",
    "design_class", "qa_flags", "duplicate_group_id", "sequence_pair_key",
    "curation_status",
]

MANUAL_REVIEW_COLUMNS = [
    "manual_target_is_antigen", "manual_label_supported_by_source",
    "manual_binder_type_classification", "manual_pairing_correct",
    "manual_notes",
]

CURATION_COLUMNS = [
    "curator", "curation_date", "curation_notes", "normalized_KD_nM",
]

AA_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

ANTIBODY_TARGET_KEYWORDS = [
    "antibody", "immunoglobulin", " fab", "fab ", "fv", "scfv", "vhh",
    "heavy chain", "light chain", "kappa chain", "lambda chain",
    "variable heavy", "variable light", "binding ig", "igm rf", "igg fc",
]


def clean(v: str | None) -> str:
    return (v or "").strip()


def clean_seq(v: str | None) -> str:
    return re.sub(r"\s+", "", clean(v).upper())


def clean_name(v: str | None) -> str:
    return normalize_antigen(clean(v))


def valid_seq(s: str) -> bool:
    return bool(s) and all(c in AA_CHARS for c in s)


def seq_key(binder_seq: str, target_seq: str) -> str:
    combined = f"{binder_seq}|||{target_seq}"
    return hashlib.md5(combined.encode()).hexdigest()


def classify_type(desc: str) -> str:
    d = desc.lower()
    d_compact = re.sub(r"[^a-z0-9]+", "", d)
    if "nanobody" in d_compact or "vhh" in d_compact:
        return "nanobody"
    if "scfv" in d_compact:
        return "scFv"
    if "heavychain" in d_compact or re.search(r"\bvh\b", d):
        return "heavy chain"
    if "lightchain" in d_compact or re.search(r"\bvl\b", d):
        return "light chain"
    if "fab" in d:
        return "Fab"
    if "cdr" in d_compact:
        return "CDR-only"
    return "other"


def merge_flags(*flag_sets: str, extra: list[str] | None = None) -> str:
    flags: list[str] = []
    for flag_set in flag_sets:
        flags.extend(flag for flag in clean(flag_set).split(";") if flag)
    flags.extend(extra or [])
    return ";".join(dict.fromkeys(flags))


def looks_like_antibody_target(desc: str) -> bool:
    d = f" {desc.lower()} "
    return any(keyword in d for keyword in ANTIBODY_TARGET_KEYWORDS)


def extract_patent_reference(row: dict) -> str:
    explicit = clean(row.get("patent_id"))
    source = clean(row.get("source"))
    text = " ".join(
        part for part in [explicit, source, clean(row.get("binder_name"))] if part
    )
    match = re.search(r"\b(?:US|WO|EP)\s*\d+[A-Z0-9 ]{0,8}", text, re.IGNORECASE)
    if match:
        return re.sub(r"\s+", " ", match.group(0)).strip(" .;,")
    if source and source.upper() != "US":
        return source
    accession = clean(row.get("ncbi_accession"))
    return f"GenBank:{accession}" if accession else source


def count_literature_rows() -> int:
    if not LITERATURE_PATH.exists():
        return 0
    with LITERATURE_PATH.open(newline="", encoding="utf-8") as f:
        return sum(1 for _ in csv.DictReader(f))


def load_seed() -> list[dict]:
    if not SEED_PATH.exists():
        print(f"  Seed not found: {SEED_PATH}")
        return []
    with SEED_PATH.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    candidates = []
    for i, row in enumerate(rows, 1):
        bseq = clean_seq(row.get("binder_sequence"))
        tseq = clean_seq(row.get("target_sequence"))
        target = clean_name(row.get("target_name"))
        btype = classify_type(f"{clean(row.get('binder_type'))} {clean(row.get('design_class'))}")
        candidates.append({
            "record_id": f"AASEED-{i:06d}",
            "binder_name": clean(row.get("binder_name")),
            "binder_type": btype,
            "binder_sequence": bseq,
            "antigen_name": target,
            "target_name": target,
            "target_sequence": tseq,
            "interaction_label": clean(row.get("interaction_label")),
            "affinity_value": clean(row.get("affinity_value")),
            "affinity_unit": clean(row.get("affinity_unit")),
            "evidence_text": clean(row.get("evidence_text")),
            "source_reference": clean(row.get("source_reference")),
            "confidence": "medium",
            "needs_review": "True",
            "source_type": clean(row.get("source_type")) or "proteinbase",
            "source_dataset": clean(row.get("source_dataset")) or "ProteinBase",
            "source_id": clean(row.get("source_id")),
            "source_confidence": clean(row.get("confidence")) or "medium",
            "design_class": clean(row.get("design_class")),
            "qa_flags": merge_flags(row.get("qa_flags", "")),
            "duplicate_group_id": "",
            "sequence_pair_key": seq_key(bseq, tseq),
            "curation_status": "candidate_unreviewed",
        })
    return candidates


def load_pdb() -> list[dict]:
    if not PDB_PATH.exists():
        print(f"  PDB data not found: {PDB_PATH}")
        return []
    with PDB_PATH.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    candidates = []
    for i, row in enumerate(rows, 1):
        bseq = clean_seq(row.get("binder_sequence"))
        aseq = clean_seq(row.get("antigen_sequence"))
        antigen = clean_name(row.get("antigen_name"))
        flags = []
        if looks_like_antibody_target(clean(row.get("antigen_name"))):
            flags.append("possible_pdb_antibody_chain_as_antigen")
        candidates.append({
            "record_id": f"AAPDB-{i:06d}",
            "binder_name": clean(row.get("binder_name")),
            "binder_type": clean(row.get("binder_type")),
            "binder_sequence": bseq,
            "antigen_name": antigen,
            "target_name": antigen,
            "target_sequence": aseq,
            "interaction_label": "binder",
            "affinity_value": "",
            "affinity_unit": "",
            "evidence_text": clean(row.get("title")),
            "source_reference": clean(row.get("source")),
            "confidence": "medium" if flags else "high",
            "needs_review": "True",
            "source_type": "PDB",
            "source_dataset": "RCSB PDB",
            "source_id": clean(row.get("pdb_id")),
            "source_confidence": "high",
            "design_class": "",
            "qa_flags": merge_flags(extra=flags),
            "duplicate_group_id": "",
            "sequence_pair_key": seq_key(bseq, aseq),
            "curation_status": "candidate_unreviewed",
        })
    return candidates


def load_patent() -> list[dict]:
    if not PATENT_PATH.exists():
        print(f"  Patent data not found: {PATENT_PATH}")
        return []
    with PATENT_PATH.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    candidates = []
    for i, row in enumerate(rows, 1):
        bseq = clean_seq(row.get("binder_sequence"))
        aseq = clean_seq(row.get("antigen_sequence"))
        antigen = clean_name(row.get("antigen_name"))
        source_ref = extract_patent_reference(row)
        candidates.append({
            "record_id": f"AAPAT-{i:06d}",
            "binder_name": clean(row.get("binder_name")),
            "binder_type": clean(row.get("binder_type")),
            "binder_sequence": bseq,
            "antigen_name": antigen,
            "target_name": antigen,
            "target_sequence": aseq,
            "interaction_label": "binder",
            "affinity_value": "",
            "affinity_unit": "",
            "evidence_text": f"NCBI Protein patent entry: {clean(row.get('binder_name'))}",
            "source_reference": source_ref,
            "confidence": "medium",
            "needs_review": "True",
            "source_type": "NCBI Patent",
            "source_dataset": "NCBI GenBank PAT",
            "source_id": clean(row.get("ncbi_accession")),
            "source_confidence": "medium",
            "design_class": "",
            "qa_flags": "",
            "duplicate_group_id": "",
            "sequence_pair_key": seq_key(bseq, aseq),
            "curation_status": "candidate_unreviewed",
        })
    return candidates


def deduplicate(rows: list[dict]) -> list[dict]:
    seen: dict[str, int] = {}
    deduped = []
    dup_count = 0
    for row in rows:
        key = row["sequence_pair_key"]
        if key in seen:
            dup_count += 1
            row["duplicate_group_id"] = f"DUP-{seen[key]:06d}"
        else:
            seen[key] = len(deduped) + 1
            deduped.append(row)
    if dup_count:
        print(f"  Removed {dup_count} duplicate sequence pairs")
    return deduped


def stratified_sample(rows: list[dict], size: int = 120, seed: int = 20260630) -> list[dict]:
    rng = random.Random(seed)
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row["source_type"], row["binder_type"])].append(row)

    per_group = max(1, size // max(1, len(groups)))
    selected: list[dict] = []
    selected_ids: set[str] = set()
    for key in sorted(groups):
        group = list(groups[key])
        rng.shuffle(group)
        for row in group[:per_group]:
            selected.append(row)
            selected_ids.add(row["record_id"])

    if len(selected) < size:
        remaining = [r for r in rows if r["record_id"] not in selected_ids]
        rng.shuffle(remaining)
        selected.extend(remaining[:size - len(selected)])

    selected = selected[:size]
    selected.sort(key=lambda r: (r["source_type"], r["binder_type"], r["antigen_name"]))

    sample = []
    for row in selected:
        entry = {col: "" for col in MANUAL_REVIEW_COLUMNS}
        entry.update(row)
        sample.append(entry)
    return sample


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_memo(candidates: list[dict], sample: list[dict]) -> None:
    source_counts = Counter(r["source_type"] for r in candidates)
    label_counts = Counter(r["interaction_label"] for r in candidates)
    type_counts = Counter(r["binder_type"] for r in candidates)
    conf_counts = Counter(r["confidence"] for r in candidates)
    antigen_counts = Counter(r["antigen_name"] for r in candidates)
    n_valid = sum(1 for r in candidates if valid_seq(r["binder_sequence"]) and valid_seq(r["target_sequence"]))
    n_affinity = sum(1 for r in candidates if r["affinity_value"])
    sample_src = Counter(r["source_type"] for r in sample)
    sample_type = Counter(r["binder_type"] for r in sample)
    sample_antigens = len(set(r["antigen_name"] for r in sample))
    top_ag = antigen_counts.most_common(8)
    literature_rows = count_literature_rows()
    generated_on = date.today().isoformat()

    text = f"""# Antibody-Antigen Validation Memo

Generated on {generated_on}. Combines sequence-pair seed data with PDB and patent extractions.

## Outputs

- `antibody_antigen_candidates.csv`: {len(candidates)} antibody-antigen candidate rows.
- `antibody_antigen_master_v2.csv`: app-facing master table ({len(candidates)} rows).
- `candidates_v1.csv`: seed-only version ({source_counts.get('ProteinBase', 0)} rows).
- `candidates_v2.csv`: all sources merged ({len(candidates)} rows).
- `antibody_antigen_curated.csv`: header-only curated table until manual review is applied.
- `manual_validation_sample.csv`: {len(sample)} rows stratified by source and binder type for manual review.
- `pubmed_references.csv`: {literature_rows} PubMed reference rows, kept as literature evidence rather than sequence-pair candidates.

## Data Sources

| Source | Rows | Confidence | Notes |
|--------|------|------------|-------|
| ProteinBase seed | {source_counts.get('ProteinBase', 0)} | medium | Sequence-positive antibody rows from protein interaction pilot |
| RCSB PDB | {source_counts.get('PDB', 0)} | high | Co-crystal structures of antibody-antigen complexes |
| NCBI Patent (GenBank) | {source_counts.get('NCBI Patent', 0)} | medium | Patent-deposited antibody sequences with annotated targets |
| **Total** | **{len(candidates)}** | | |

PubMed is not counted as a sequence-pair source because abstracts usually do not
provide paired binder and antigen amino-acid sequences. It is used as a reference
and evidence layer.

## Candidate QA

- Total rows: {len(candidates)}
- Rows with valid binder and target sequences: {n_valid}/{len(candidates)}
- Rows with affinity values: {n_affinity}/{len(candidates)}
- Interaction label: {', '.join(f'{k} {v}' for k, v in label_counts.most_common())}
- Binder types: {', '.join(f'{k} {v}' for k, v in type_counts.most_common())}
- Confidence: {', '.join(f'{k} {v}' for k, v in conf_counts.most_common())}
- Unique antigens: {len(antigen_counts)}

### Top Antigens

{chr(10).join(f'- {name}: {count} ({count * 100 // len(candidates)}%)' for name, count in top_ag)}

## Manual Sample

{len(sample)} rows stratified by source x binder type:
- Sources: {', '.join(f'{k}: {v}' for k, v in sample_src.most_common())}
- Types: {', '.join(f'{k}: {v}' for k, v in sample_type.most_common())}
- Covers {sample_antigens} unique antigens

## Pipeline

```
scripts/fetch_pdb_complexes.py   -> data/pdb/pdb_ab_ag_pairs.csv
scripts/fetch_patent_sequences.py -> data/patent/patent_ab_ag_pairs.csv
scripts/fetch_pubmed_references.py -> pubmed_references.csv
scripts/normalize_antigens.py     (normalizes antigen names in-place)
scripts/merge_candidates.py       -> antibody_antigen_master_v2.csv
                                  -> outputs/candidates_v1.csv, candidates_v2.csv
scripts/build_antibody_antigen_candidates.py  (original seed-only builder)
app.py                            (Streamlit browser)
```

## Next Steps

1. Reclassify remaining "other" + "antibody fragment" rows during manual review.
2. Further normalize PDB-derived antigen names to standard gene symbols.
3. Add published KD values from literature for PDB-derived pairs.
4. Expand NCBI searches to additional targets (IL-17, IL-4/IL-13, PCSK9, CD38, CGRP, BCMA, TROP2, etc.).
5. Promote reviewed rows into `antibody_antigen_curated.csv`.
"""
    memo_path = OUTPUT_DIR / "validation_memo.md"
    memo_path.write_text(text, encoding="utf-8")
    print(f"  Wrote {memo_path}")


def main() -> None:
    print("Loading sources...")
    seed = load_seed()
    pdb = load_pdb()
    patent = load_patent()
    print(f"  Seed: {len(seed)}, PDB: {len(pdb)}, Patent: {len(patent)}")

    all_rows = seed + pdb + patent
    print(f"Total before dedup: {len(all_rows)}")
    all_rows = deduplicate(all_rows)
    print(f"Total after dedup: {len(all_rows)}")

    write_csv(OUTPUT_DIR / "candidates_v1.csv", CANDIDATE_COLUMNS, seed)
    print(f"  candidates_v1.csv: {len(seed)} rows (seed only)")

    write_csv(OUTPUT_DIR / "candidates_v2.csv", CANDIDATE_COLUMNS, all_rows)
    print(f"  candidates_v2.csv: {len(all_rows)} rows (all sources)")

    write_csv(OUTPUT_DIR / "antibody_antigen_candidates.csv", CANDIDATE_COLUMNS, all_rows)
    print(f"  antibody_antigen_candidates.csv: {len(all_rows)} rows")

    write_csv(MASTER_PATH, CANDIDATE_COLUMNS, all_rows)
    print(f"  antibody_antigen_master_v2.csv: {len(all_rows)} rows")

    curated_cols = CANDIDATE_COLUMNS[:14] + MANUAL_REVIEW_COLUMNS + CURATION_COLUMNS + CANDIDATE_COLUMNS[14:]
    write_csv(OUTPUT_DIR / "antibody_antigen_curated.csv", curated_cols, [])

    sample = stratified_sample(all_rows)
    write_csv(OUTPUT_DIR / "manual_validation_sample.csv", MANUAL_REVIEW_COLUMNS + CANDIDATE_COLUMNS, sample)
    print(f"  manual_validation_sample.csv: {len(sample)} rows")

    write_memo(all_rows, sample)
    print("\nDone.")


if __name__ == "__main__":
    main()
