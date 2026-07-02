#!/usr/bin/env python3
"""
Build the standalone antibody-antigen candidate files.

This project starts from the sequence-positive antibody rows copied out of the
protein interaction pilot. It intentionally does not mix in the de novo binder
track; that file is kept only as a reference under data/source_pilots/.
"""

from __future__ import annotations

import csv
import random
import re
from collections import Counter, defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PILOT = PROJECT_ROOT / "data" / "source_pilots" / "antibody_antigen_pilot.csv"
OUTPUT_DIR = PROJECT_ROOT / "outputs"

CANDIDATES_PATH = OUTPUT_DIR / "antibody_antigen_candidates.csv"
CURATED_PATH = OUTPUT_DIR / "antibody_antigen_curated.csv"
MANUAL_SAMPLE_PATH = OUTPUT_DIR / "manual_validation_sample.csv"
MEMO_PATH = OUTPUT_DIR / "validation_memo.md"

AA_EXTENDED = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

CORE_COLUMNS = [
    "record_id",
    "binder_name",
    "binder_type",
    "binder_sequence",
    "antigen_name",
    "target_name",
    "target_sequence",
    "interaction_label",
    "affinity_value",
    "affinity_unit",
    "evidence_text",
    "source_reference",
    "confidence",
    "needs_review",
]

PROVENANCE_COLUMNS = [
    "source_type",
    "source_dataset",
    "source_id",
    "source_confidence",
    "design_class",
    "qa_flags",
    "duplicate_group_id",
    "sequence_pair_key",
    "curation_status",
]

CANDIDATE_COLUMNS = CORE_COLUMNS + PROVENANCE_COLUMNS

MANUAL_REVIEW_COLUMNS = [
    "manual_target_is_antigen",
    "manual_label_supported_by_source",
    "manual_binder_type_classification",
    "manual_pairing_correct",
    "manual_notes",
]

CURATION_COLUMNS = [
    "curator",
    "curation_date",
    "curation_notes",
    "normalized_KD_nM",
]

CURATED_COLUMNS = CORE_COLUMNS + MANUAL_REVIEW_COLUMNS + CURATION_COLUMNS + PROVENANCE_COLUMNS


def clean_text(value: str | None) -> str:
    return (value or "").strip()


def clean_sequence(value: str | None) -> str:
    return re.sub(r"\s+", "", clean_text(value).upper())


def valid_sequence(sequence: str) -> bool:
    return bool(sequence) and all(char in AA_EXTENDED for char in sequence)


def normalize_name(value: str | None) -> str:
    text = clean_text(value)
    text = re.sub(r"\s+", " ", text)
    return text


def standardize_binder_type(value: str | None, design_class: str | None = "") -> str:
    text = f"{clean_text(value)} {clean_text(design_class)}".lower()
    compact = re.sub(r"[^a-z0-9]+", "", text)
    if "nanobody" in compact or "vhh" in compact:
        return "nanobody"
    if "scfv" in compact:
        return "scFv"
    if "heavychain" in compact or re.search(r"\bvh\b", text):
        return "heavy chain"
    if "lightchain" in compact or re.search(r"\bvl\b", text):
        return "light chain"
    if "cdr" in compact:
        return "CDR-only"
    return "other"


def merge_flags(*flag_sets: str, extra: list[str] | None = None) -> str:
    flags: list[str] = []
    for flag_set in flag_sets:
        flags.extend(flag for flag in clean_text(flag_set).split(";") if flag)
    flags.extend(extra or [])
    return ";".join(dict.fromkeys(flags))


def load_source_rows() -> list[dict[str, str]]:
    with SOURCE_PILOT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_candidates(source_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for index, row in enumerate(source_rows, start=1):
        binder_sequence = clean_sequence(row.get("binder_sequence"))
        target_sequence = clean_sequence(row.get("target_sequence"))
        binder_type = standardize_binder_type(row.get("binder_type"), row.get("design_class"))

        project_flags = []
        if not valid_sequence(binder_sequence):
            project_flags.append("project_invalid_binder_sequence")
        if not valid_sequence(target_sequence):
            project_flags.append("project_invalid_target_sequence")
        if binder_type == "other":
            project_flags.append("project_antibody_fragment_type_unresolved")
        project_flags.append("project_antigen_role_unreviewed")
        project_flags.append("project_label_support_unreviewed")

        candidate = {
            "record_id": f"AASEED-{index:06d}",
            "binder_name": normalize_name(row.get("binder_name")),
            "binder_type": binder_type,
            "binder_sequence": binder_sequence,
            "antigen_name": normalize_name(row.get("target_name")),
            "target_name": normalize_name(row.get("target_name")),
            "target_sequence": target_sequence,
            "interaction_label": clean_text(row.get("interaction_label")),
            "affinity_value": clean_text(row.get("affinity_value")),
            "affinity_unit": clean_text(row.get("affinity_unit")),
            "evidence_text": clean_text(row.get("evidence_text")),
            "source_reference": clean_text(row.get("source_reference")),
            "confidence": "medium",
            "needs_review": "True",
            "source_type": clean_text(row.get("source_type")),
            "source_dataset": clean_text(row.get("source_dataset")),
            "source_id": clean_text(row.get("source_id")),
            "source_confidence": clean_text(row.get("confidence")),
            "design_class": clean_text(row.get("design_class")),
            "qa_flags": merge_flags(row.get("qa_flags", ""), extra=project_flags),
            "duplicate_group_id": clean_text(row.get("duplicate_group_id")),
            "sequence_pair_key": clean_text(row.get("sequence_pair_key")),
            "curation_status": "candidate_unreviewed",
        }
        candidates.append(candidate)
    return candidates


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def stratified_manual_sample(rows: list[dict[str, str]], size: int = 50, seed: int = 20260629) -> list[dict[str, str]]:
    if len(rows) <= size:
        selected = list(rows)
    else:
        rng = random.Random(seed)
        groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            groups[(row["binder_type"], row["interaction_label"])].append(row)

        selected_ids: set[str] = set()
        selected = []
        per_group = max(1, size // max(1, len(groups)))
        for key in sorted(groups):
            group = list(groups[key])
            rng.shuffle(group)
            for row in group[: min(per_group, len(group))]:
                selected.append(row)
                selected_ids.add(row["record_id"])

        if len(selected) < size:
            remaining = [row for row in rows if row["record_id"] not in selected_ids]
            rng.shuffle(remaining)
            selected.extend(remaining[: size - len(selected)])

        selected = selected[:size]

    selected.sort(key=lambda row: (row["binder_type"], row["interaction_label"], row["antigen_name"], row["binder_name"]))
    sample = []
    for row in selected:
        review_fields = {column: "" for column in MANUAL_REVIEW_COLUMNS}
        sample.append(review_fields | row)
    return sample


def count_nonempty(rows: list[dict[str, str]], column: str) -> int:
    return sum(1 for row in rows if clean_text(row.get(column)))


def count_valid_pairs(rows: list[dict[str, str]]) -> int:
    return sum(
        1
        for row in rows
        if valid_sequence(row.get("binder_sequence", "")) and valid_sequence(row.get("target_sequence", ""))
    )


def format_counts(counter: Counter[str]) -> str:
    if not counter:
        return "none"
    return ", ".join(f"{key}: {counter[key]}" for key in sorted(counter))


def compute_antigen_stats(candidates: list[dict[str, str]]) -> tuple[int, Counter[str]]:
    antigen_same = sum(1 for r in candidates if r["antigen_name"] == r["target_name"])
    antigen_counts = Counter(r["antigen_name"] for r in candidates)
    return antigen_same, antigen_counts


def compute_affinity_stats(candidates: list[dict[str, str]]) -> tuple[int, Counter[str], str]:
    unit_counts = Counter(r["affinity_unit"] for r in candidates if r["affinity_unit"])
    vals = [float(r["affinity_value"]) for r in candidates if r["affinity_value"]]
    if vals:
        range_str = f"{min(vals):.2e} to {max(vals):.2e}"
    else:
        range_str = "n/a"
    return len(vals), unit_counts, range_str


def write_memo(candidates: list[dict[str, str]], manual_sample: list[dict[str, str]]) -> None:
    source_counts = Counter(row["source_type"] for row in candidates)
    label_counts = Counter(row["interaction_label"] for row in candidates)
    binder_type_counts = Counter(row["binder_type"] for row in candidates)
    source_confidence_counts = Counter(row["source_confidence"] for row in candidates)
    project_confidence_counts = Counter(row["confidence"] for row in candidates)
    needs_review_counts = Counter(row["needs_review"] for row in candidates)

    antigen_same, antigen_counts = compute_antigen_stats(candidates)
    n_affinity, affinity_units, affinity_range = compute_affinity_stats(candidates)
    top2 = antigen_counts.most_common(2)
    top2_pct = sum(v for _, v in top2) * 100 // len(candidates) if candidates else 0
    top2_str = ", ".join(f"{name} {count}/{len(candidates)} ({count * 100 // len(candidates)}%)" for name, count in top2)

    manual_antigen_counts = Counter(r["antigen_name"] for r in manual_sample)
    manual_coverage = f"{len(manual_antigen_counts)}/{len(antigen_counts)}"

    text = f"""# Antibody-Antigen Validation Memo

Generated by `scripts/build_antibody_antigen_candidates.py`.

## Outputs

- `antibody_antigen_candidates.csv`: {len(candidates)} sequence-positive antibody-antigen candidate rows.
- `antibody_antigen_curated.csv`: header-only curated table until manual review is applied.
- `manual_validation_sample.csv`: {len(manual_sample)} antibody-only rows for manual review.

## Input Provenance

- Source seed: `data/source_pilots/antibody_antigen_pilot.csv`
- De novo reference kept separate: `data/source_pilots/de_novo_binder_pilot_reference.csv`
- Raw copied sources: `data/raw/proteinbase_interactions.csv`, `data/raw/litscrape.csv`
- Original project copied from `/Users/chen/reprogramming_project`; no source files were moved.

## Candidate QA

- Rows with valid binder and target amino-acid sequences: {count_valid_pairs(candidates)}/{len(candidates)}
- Rows with affinity values: {n_affinity}/{len(candidates)} (all tagged {format_counts(affinity_units)}; range {affinity_range}, likely in molar units)
- Source counts: {format_counts(source_counts)}
- Interaction label counts: {format_counts(label_counts)}
- Binder type counts: {format_counts(binder_type_counts)}
- Source confidence counts: {format_counts(source_confidence_counts)}
- Project confidence counts: {format_counts(project_confidence_counts)}
- Project needs_review counts: {format_counts(needs_review_counts)}

### Antigen-Name Status

`antigen_name` is currently an exact copy of `target_name` for {antigen_same}/{len(candidates)} rows. The
two fields exist to allow correction during manual review — some targets may be
receptors, enzymes, or other proteins where the immunological antigen is a different
entity (e.g., the target is a receptor but the antigen used for immunization was a
soluble ectodomain or a peptide fragment). Reviewers should verify whether
`antigen_name` needs to diverge from `target_name` and update accordingly.

### Antigen Concentration

Top 2 antigens account for {top2_pct}% of candidates: {top2_str}. This reflects
the source seed composition. Patent-scale extraction should diversify the antigen
distribution; until then, per-antigen metrics should be reported alongside aggregate
numbers to avoid masking low-coverage antigens.

Manual sample covers {manual_coverage} unique antigens from the candidate set.

## Interpretation

This is an antibody-antigen seed, not a finished patent/literature extraction.
The rows have paired binder and target amino-acid sequences, but the project marks
every row as `needs_review=True` until manual review confirms that the target is
an antigen, the binder-target label is supported by the source, and the antibody
fragment type is correct.

## Manual Review Fields

- `manual_target_is_antigen`: confirm that `target_name`/`antigen_name` is really the antigen.
- `manual_label_supported_by_source`: confirm binder/non-binder label against the original source.
- `manual_binder_type_classification`: choose nanobody, scFv, heavy chain, light chain, CDR-only, or other.
- `manual_pairing_correct`: confirm sequence-to-antigen pairing.
- `manual_notes`: capture source-specific concerns or evidence snippets.

## Patent-Scale Extraction Next Steps

1. Parse patent sequence listings and normalize every `SEQ ID NO`.
2. Classify antibody fragments as heavy chain, light chain, paired heavy/light, scFv, nanobody, CDR-only, or other.
3. Map each sequence to antigen names and evidence sentences.
4. Resolve heavy/light-chain pairing from claims, examples, tables, and sequence-listing context.
5. Promote rows into `antibody_antigen_curated.csv` only after source evidence and pairing are checked.
"""
    MEMO_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    source_rows = load_source_rows()
    candidates = build_candidates(source_rows)
    manual_sample = stratified_manual_sample(candidates)

    write_csv(CANDIDATES_PATH, CANDIDATE_COLUMNS, candidates)
    write_csv(CURATED_PATH, CURATED_COLUMNS, [])
    write_csv(MANUAL_SAMPLE_PATH, MANUAL_REVIEW_COLUMNS + CANDIDATE_COLUMNS, manual_sample)
    write_memo(candidates, manual_sample)

    print(f"Wrote {CANDIDATES_PATH} ({len(candidates)} rows)")
    print(f"Wrote {CURATED_PATH} (0 curated rows)")
    print(f"Wrote {MANUAL_SAMPLE_PATH} ({len(manual_sample)} rows)")
    print(f"Wrote {MEMO_PATH}")


if __name__ == "__main__":
    main()
