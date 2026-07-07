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
    "antigen_uniprot_id", "binder_species", "antigen_species",
    "interaction_label", "affinity_value", "affinity_unit", "normalized_KD_nM",
    "evidence_text", "source_reference", "confidence", "needs_review",
    "source_type", "source_dataset", "source_id", "source_confidence",
    "design_class", "complex_group_id", "pdb_pairing_status",
    "paired_chain_roles", "binder_type_status", "binder_type_suggestion",
    "qa_flags", "duplicate_group_id", "sequence_pair_key", "curation_status",
]

MANUAL_REVIEW_COLUMNS = [
    "manual_target_is_antigen", "manual_label_supported_by_source",
    "manual_binder_type_classification", "manual_pairing_correct",
    "manual_notes",
]

CURATION_COLUMNS = [
    "curator", "curation_date", "curation_notes",
]

AA_STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
AA_AMBIGUOUS = set("BXZUO")
AA_CHARS = AA_STANDARD | AA_AMBIGUOUS

PATENT_ANTIGEN_UNIPROT = {
    "PD-1": "Q15116", "PD-L1": "Q9NZQ7", "CTLA-4": "P16410",
    "TIGIT": "Q495A1", "LAG-3": "P18627", "TIM-3": "Q8TDQ0",
    "VISTA": "Q9H7M9", "BTLA": "Q7Z6A9", "CD47": "Q08722", "SIRPa": "P78324",
    "HER2": "P04626", "EGFR": "P00533", "HER3": "P21860",
    "CD20": "P11836", "CD19": "P15391", "CD3": "P07766",
    "CD38": "P28907", "BCMA": "Q02223", "GPRC5D": "Q9NZD1", "FcRH5": "Q96RD9",
    "CD33": "P20138", "CD22": "P20273", "CD30": "P28908",
    "TROP2": "P09758", "Nectin-4": "Q96NY8", "CLDN18.2": "P56856",
    "DLL3": "Q9NYJ7", "Mesothelin": "Q13421", "MET": "P08581",
    "PSMA": "Q04609", "FAP": "Q12884", "GPC3": "P51654",
    "B7-H3": "Q5ZPR3", "EpCAM": "P16422",
    "TNF-alpha": "P01375", "IL-6": "P05231", "IL-6R": "P08887",
    "IL-17A": "Q16552", "IL-17RA": "Q96F46", "IL-4Ra": "P24394",
    "IL-13": "P35225", "IL-5": "P05113", "IL-5Ra": "Q01344",
    "IL-23": "Q9NPF7", "IL-12/23 p40": "P29460",
    "IL-33": "O95760", "IL-31": "Q6EBC2", "TSLP": "Q969D9",
    "IgE": "P01854", "VEGF": "P15692", "VEGFR2": "P35968",
    "NGF": "P01138", "CGRP": "P06881", "CGRPr": "Q16602",
    "PCSK9": "Q8NBP7", "C5": "P01031", "RANKL": "O14788",
    "Ang2": "O15123", "BAFF": "Q9Y275",
    "RSV F": "P03420", "RSV": "P03420",
    "SARS-CoV-2 Spike": "P0DTC2",
    "Influenza HA": "P03437", "HIV gp120": "P04578", "Ebola GP": "Q05320",
}

ANTIGEN_SPECIES_BY_NAME = {
    "PD-1": "Homo sapiens", "PD-L1": "Homo sapiens", "CTLA-4": "Homo sapiens",
    "TIGIT": "Homo sapiens", "LAG-3": "Homo sapiens", "TIM-3": "Homo sapiens",
    "VISTA": "Homo sapiens", "BTLA": "Homo sapiens", "CD47": "Homo sapiens",
    "SIRPa": "Homo sapiens",
    "HER2": "Homo sapiens", "EGFR": "Homo sapiens", "HER3": "Homo sapiens",
    "CD20": "Homo sapiens", "CD19": "Homo sapiens", "CD3": "Homo sapiens",
    "CD38": "Homo sapiens", "BCMA": "Homo sapiens", "GPRC5D": "Homo sapiens",
    "FcRH5": "Homo sapiens", "CD33": "Homo sapiens", "CD22": "Homo sapiens",
    "CD30": "Homo sapiens", "TROP2": "Homo sapiens", "Nectin-4": "Homo sapiens",
    "CLDN18.2": "Homo sapiens", "DLL3": "Homo sapiens", "Mesothelin": "Homo sapiens",
    "MET": "Homo sapiens", "PSMA": "Homo sapiens", "FAP": "Homo sapiens",
    "GPC3": "Homo sapiens", "B7-H3": "Homo sapiens", "EpCAM": "Homo sapiens",
    "TNF-alpha": "Homo sapiens", "IL-6": "Homo sapiens", "IL-6R": "Homo sapiens",
    "IL-17A": "Homo sapiens", "IL-17RA": "Homo sapiens", "IL-4Ra": "Homo sapiens",
    "IL-13": "Homo sapiens", "IL-5": "Homo sapiens", "IL-5Ra": "Homo sapiens",
    "IL-23": "Homo sapiens", "IL-12/23 p40": "Homo sapiens",
    "IL-33": "Homo sapiens", "IL-31": "Homo sapiens", "TSLP": "Homo sapiens",
    "IgE": "Homo sapiens", "VEGF": "Homo sapiens", "VEGFR2": "Homo sapiens",
    "NGF": "Homo sapiens", "CGRP": "Homo sapiens", "CGRPr": "Homo sapiens",
    "PCSK9": "Homo sapiens", "C5": "Homo sapiens", "RANKL": "Homo sapiens",
    "Ang2": "Homo sapiens", "BAFF": "Homo sapiens",
    "Human serum albumin": "Homo sapiens",
    "IL-7R": "Homo sapiens", "RBX1": "Homo sapiens",
    "hnmt": "Homo sapiens", "human ambp": "Homo sapiens",
    "human gm2a": "Homo sapiens", "human idi2": "Homo sapiens",
    "human insulin receptor": "Homo sapiens",
    "human mzb1 perp1": "Homo sapiens", "human orm2": "Homo sapiens",
    "human pdgfr beta": "Homo sapiens", "human phyh": "Homo sapiens",
    "human pmvk": "Homo sapiens", "human rfk": "Homo sapiens",
    "RSV F": "Respiratory syncytial virus",
    "RSV": "Respiratory syncytial virus",
    "SARS-CoV-2 Spike": "Severe acute respiratory syndrome coronavirus 2",
    "Influenza HA": "Influenza A virus",
    "HIV gp120": "Human immunodeficiency virus 1",
    "Ebola GP": "Zaire ebolavirus",
    "Nipah glycoprotein G": "Nipah virus",
}

NON_ANTIBODY_BINDER_PATTERNS = [
    r"\bprotein a\b",
    r"\bprotein l\b",
    r"\bimmunoglobulin g binding protein\b",
    r"\bimmunoglobulin-binding protein\b",
]

PATENT_NUMBER_RE = re.compile(
    r"\b(?:US|WO|EP)\s*\d{4,12}\s*(?:[A-Z]\d?|[A-Z]{1,2}\d?)?\b",
    re.IGNORECASE,
)

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


def sequence_quality_flags(seq: str, label: str) -> list[str]:
    flags: list[str] = []
    if not seq:
        flags.append(f"{label}_sequence_missing")
        return flags
    invalid = sorted({c for c in seq if c not in AA_CHARS})
    if invalid:
        flags.append(f"{label}_sequence_invalid_chars")
    if "X" in seq:
        flags.append(f"{label}_sequence_contains_x")
    ambiguous = sorted({c for c in seq if c in AA_AMBIGUOUS and c != "X"})
    if ambiguous:
        flags.append(f"{label}_sequence_ambiguous_aa")
    return flags


def row_quality_flags(binder_seq: str, target_seq: str) -> list[str]:
    return (
        sequence_quality_flags(binder_seq, "binder")
        + sequence_quality_flags(target_seq, "target")
    )


def normalized_kd_nm(value: str | None, unit: str | None) -> str:
    raw_value = clean(value)
    raw_unit = clean(unit).lower()
    if not raw_value:
        return ""
    try:
        number = float(raw_value)
    except ValueError:
        return ""

    # ProteinBase stores raw KD as molar values and labels the unit as source_KD.
    if raw_unit in {"source_kd", "m", "molar", "mol/l", "mol/liter"}:
        return f"{number * 1e9:.6g}"
    if raw_unit in {"nm", "nanomolar", "nanomole", "nanomoles"}:
        return f"{number:.6g}"
    if raw_unit in {"um", "µm", "micromolar"}:
        return f"{number * 1e3:.6g}"
    if raw_unit in {"pm", "picomolar"}:
        return f"{number * 1e-3:.6g}"
    return ""


def seq_key(binder_seq: str, target_seq: str) -> str:
    combined = f"{binder_seq}|||{target_seq}"
    return hashlib.md5(combined.encode()).hexdigest()


def classify_type(desc: str) -> str:
    d = desc.lower()
    d_compact = re.sub(r"[^a-z0-9]+", "", d)
    if (
        "nanobody" in d_compact
        or "vhh" in d_compact
        or "singledomain" in d_compact
        or "sdab" in d_compact
        or "vnar" in d_compact
        or "newantigenreceptor" in d_compact
    ):
        return "nanobody"
    if "scfv" in d_compact or "singlechainfv" in d_compact:
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


def classify_with_fallback(raw_type: str, *text_parts: str) -> str:
    text = " ".join([raw_type, *text_parts])
    inferred = classify_type(text)
    if inferred != "other":
        return inferred
    raw_clean = clean(raw_type)
    return raw_clean if raw_clean else "other"


def is_non_antibody_binder(text: str) -> bool:
    d = text.lower()
    return any(re.search(pattern, d) for pattern in NON_ANTIBODY_BINDER_PATTERNS)


def normalize_patent_number(value: str) -> str:
    match = PATENT_NUMBER_RE.search(value or "")
    if not match:
        return clean(value)
    raw = re.sub(r"\s+", " ", match.group(0).upper()).strip(" .;,")
    return raw


def infer_antigen_species(antigen_name: str, evidence_text: str = "") -> str:
    antigen = clean(antigen_name)
    if antigen in ANTIGEN_SPECIES_BY_NAME:
        return ANTIGEN_SPECIES_BY_NAME[antigen]
    text = f"{antigen} {evidence_text}".lower()
    if "sars-cov-2" in text or "sars cov 2" in text:
        return "Severe acute respiratory syndrome coronavirus 2"
    if "respiratory syncytial" in text or re.search(r"\brsv\b", text):
        return "Respiratory syncytial virus"
    if "hiv-1" in text or re.search(r"\bhiv\b", text):
        return "Human immunodeficiency virus 1"
    if "influenza" in text:
        return "Influenza virus"
    if "nipah" in text:
        return "Nipah virus"
    if "hen egg" in text or "chicken" in text:
        return "Gallus gallus"
    if "turkey" in text:
        return "Meleagris gallopavo"
    if "guinea fowl" in text:
        return "Numida meleagris"
    if "human" in text:
        return "Homo sapiens"
    if "mouse" in text or "murine" in text:
        return "Mus musculus"
    if "rat" in text:
        return "Rattus norvegicus"
    return ""


def infer_binder_species(binder_name: str, evidence_text: str = "") -> str:
    text = f"{binder_name} {evidence_text}".lower()
    if "new antigen receptor" in text or "vnar" in text:
        return "Chondrichthyes"
    if "camelid" in text:
        return "Camelidae"
    if "dromedary" in text:
        return "Camelus dromedarius"
    if "llama" in text:
        return "Lama glama"
    if "human ig" in text or "human antibody" in text or "human immunoglobulin" in text:
        return "Homo sapiens"
    if "mouse antibody" in text or "murine antibody" in text:
        return "Mus musculus"
    return ""


def binder_type_status(
    binder_type: str,
    binder_name: str,
    binder_sequence: str,
    evidence_text: str,
) -> tuple[str, str, list[str]]:
    text = f"{binder_type} {binder_name}"
    if is_non_antibody_binder(text):
        return (
            "non_antibody_candidate",
            "Review or exclude: binder is an immunoglobulin-binding protein, not an antibody domain.",
            ["non_antibody_binder_candidate"],
        )

    if binder_type in {"nanobody", "scFv", "heavy chain", "light chain", "Fab", "CDR-only"}:
        return ("classified", "", [])

    seq_len = len(binder_sequence)
    if 80 <= seq_len <= 140:
        suggestion = "Variable-domain-sized sequence; run ANARCI/manual review to separate VH, VL, or VHH."
    elif 180 <= seq_len <= 280:
        suggestion = "Fragment-sized sequence; inspect source for Fab/scFv/domain boundaries."
    elif seq_len > 280:
        suggestion = "Long antibody-fragment sequence; inspect patent SEQ ID for chain pairing or concatenation."
    else:
        suggestion = "Unresolved antibody type; inspect source metadata and sequence annotation."
    return ("needs_manual_review", suggestion, ["binder_type_unresolved"])


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
    if explicit:
        return normalize_patent_number(explicit)
    source = clean(row.get("source"))
    text = " ".join(part for part in [source, clean(row.get("binder_name"))] if part)
    match = PATENT_NUMBER_RE.search(text)
    if match:
        return normalize_patent_number(match.group(0))
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
        affinity_value = clean(row.get("affinity_value"))
        affinity_unit = clean(row.get("affinity_unit"))
        evidence_text = clean(row.get("evidence_text"))
        binder_name = clean(row.get("binder_name"))
        type_status, type_suggestion, type_flags = binder_type_status(
            btype, binder_name, bseq, evidence_text
        )
        candidates.append({
            "record_id": f"AASEED-{i:06d}",
            "binder_name": binder_name,
            "binder_type": btype,
            "binder_sequence": bseq,
            "antigen_name": target,
            "target_name": target,
            "target_sequence": tseq,
            "antigen_uniprot_id": "",
            "binder_species": "",
            "antigen_species": infer_antigen_species(target, evidence_text),
            "interaction_label": clean(row.get("interaction_label")),
            "affinity_value": affinity_value,
            "affinity_unit": affinity_unit,
            "normalized_KD_nM": normalized_kd_nm(affinity_value, affinity_unit),
            "evidence_text": evidence_text,
            "source_reference": clean(row.get("source_reference")),
            "confidence": "medium",
            "needs_review": "True",
            "source_type": clean(row.get("source_type")) or "proteinbase",
            "source_dataset": clean(row.get("source_dataset")) or "ProteinBase",
            "source_id": clean(row.get("source_id")),
            "source_confidence": clean(row.get("confidence")) or "medium",
            "design_class": clean(row.get("design_class")),
            "complex_group_id": "",
            "pdb_pairing_status": "",
            "paired_chain_roles": "",
            "binder_type_status": type_status,
            "binder_type_suggestion": type_suggestion,
            "qa_flags": merge_flags(
                row.get("qa_flags", ""),
                extra=type_flags + row_quality_flags(bseq, tseq),
            ),
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
        evidence_text = clean(row.get("title"))
        binder_name = clean(row.get("binder_name"))
        if is_non_antibody_binder(binder_name):
            btype = clean(row.get("binder_type")) or "other"
        else:
            btype = classify_with_fallback(clean(row.get("binder_type")), binder_name, evidence_text)
        flags = []
        if looks_like_antibody_target(clean(row.get("antigen_name"))):
            flags.append("possible_pdb_antibody_chain_as_antigen")
        type_status, type_suggestion, type_flags = binder_type_status(
            btype, binder_name, bseq, evidence_text
        )
        candidates.append({
            "record_id": f"AAPDB-{i:06d}",
            "binder_name": binder_name,
            "binder_type": btype,
            "binder_sequence": bseq,
            "antigen_name": antigen,
            "target_name": antigen,
            "target_sequence": aseq,
            "antigen_uniprot_id": "",
            "binder_species": infer_binder_species(binder_name, evidence_text),
            "antigen_species": infer_antigen_species(antigen, evidence_text),
            "interaction_label": "binder",
            "affinity_value": "",
            "affinity_unit": "",
            "normalized_KD_nM": "",
            "evidence_text": evidence_text,
            "source_reference": clean(row.get("source")),
            "confidence": "medium" if flags else "high",
            "needs_review": "True",
            "source_type": "PDB",
            "source_dataset": "RCSB PDB",
            "source_id": clean(row.get("pdb_id")),
            "source_confidence": "high",
            "design_class": "",
            "complex_group_id": "",
            "pdb_pairing_status": "",
            "paired_chain_roles": "",
            "binder_type_status": type_status,
            "binder_type_suggestion": type_suggestion,
            "qa_flags": merge_flags(extra=flags + type_flags + row_quality_flags(bseq, aseq)),
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
        antigen_uniprot_id = clean(row.get("antigen_uniprot_id")) or PATENT_ANTIGEN_UNIPROT.get(antigen, "")
        source_ref = extract_patent_reference(row)
        binder_name = clean(row.get("binder_name"))
        evidence_text = f"NCBI Protein patent entry: {binder_name}"
        btype = classify_with_fallback(clean(row.get("binder_type")), binder_name)
        type_status, type_suggestion, type_flags = binder_type_status(
            btype, binder_name, bseq, evidence_text
        )
        candidates.append({
            "record_id": f"AAPAT-{i:06d}",
            "binder_name": binder_name,
            "binder_type": btype,
            "binder_sequence": bseq,
            "antigen_name": antigen,
            "target_name": antigen,
            "target_sequence": aseq,
            "antigen_uniprot_id": antigen_uniprot_id,
            "binder_species": "",
            "antigen_species": infer_antigen_species(antigen),
            "interaction_label": "binder",
            "affinity_value": "",
            "affinity_unit": "",
            "normalized_KD_nM": "",
            "evidence_text": evidence_text,
            "source_reference": source_ref,
            "confidence": "medium",
            "needs_review": "True",
            "source_type": "NCBI Patent",
            "source_dataset": "NCBI GenBank PAT",
            "source_id": clean(row.get("ncbi_accession")),
            "source_confidence": "medium",
            "design_class": "",
            "complex_group_id": "",
            "pdb_pairing_status": "",
            "paired_chain_roles": "",
            "binder_type_status": type_status,
            "binder_type_suggestion": type_suggestion,
            "qa_flags": merge_flags(extra=type_flags + row_quality_flags(bseq, aseq)),
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


def role_summary(rows: list[dict]) -> str:
    roles = Counter(r["binder_type"] for r in rows if r.get("binder_type"))
    return ";".join(f"{role}:{count}" for role, count in sorted(roles.items()))


def pdb_pairing_status(row: dict, group_rows: list[dict]) -> str:
    roles = Counter(r["binder_type"] for r in group_rows)
    btype = row.get("binder_type", "")
    has_heavy = roles.get("heavy chain", 0) > 0
    has_light = roles.get("light chain", 0) > 0
    if has_heavy and has_light and btype in {"heavy chain", "light chain"}:
        return "vh_vl_pair_available"
    if btype == "heavy chain":
        return "heavy_chain_without_light_in_group"
    if btype == "light chain":
        return "light_chain_without_heavy_in_group"
    if btype == "scFv":
        return "single_chain_scfv"
    if btype == "nanobody":
        return "single_domain_binder"
    if btype == "Fab":
        return "fab_or_fab_like_fragment"
    if row.get("binder_type_status") == "non_antibody_candidate":
        return "non_antibody_binder_candidate"
    return "pdb_pairing_unresolved"


def annotate_pdb_pairing(rows: list[dict]) -> None:
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("source_type") == "PDB":
            groups[(row.get("source_id", ""), row.get("antigen_name", ""))].append(row)

    for (pdb_id, antigen), group_rows in groups.items():
        summary = role_summary(group_rows)
        group_id = f"PDB:{pdb_id}|{antigen}"
        for row in group_rows:
            row["complex_group_id"] = group_id
            row["paired_chain_roles"] = summary
            row["pdb_pairing_status"] = pdb_pairing_status(row, group_rows)


REVIEW_QUEUE_COLUMNS = [
    "record_id", "source_type", "source_id", "binder_name", "binder_type",
    "binder_type_status", "binder_type_suggestion", "binder_sequence_length",
    "antigen_name", "source_reference", "qa_flags",
]


def binder_type_review_queue(rows: list[dict]) -> list[dict]:
    queue = []
    for row in rows:
        if row.get("binder_type_status") == "classified":
            continue
        entry = {col: row.get(col, "") for col in REVIEW_QUEUE_COLUMNS}
        entry["binder_sequence_length"] = str(len(row.get("binder_sequence", "")))
        queue.append(entry)
    queue.sort(key=lambda r: (r["binder_type_status"], r["source_type"], r["binder_type"], r["record_id"]))
    return queue


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
    n_normalized_kd = sum(1 for r in candidates if r.get("normalized_KD_nM"))
    n_uniprot = sum(1 for r in candidates if r.get("antigen_uniprot_id"))
    n_antigen_species = sum(1 for r in candidates if r.get("antigen_species"))
    n_binder_species = sum(1 for r in candidates if r.get("binder_species"))
    pairing_counts = Counter(r.get("pdb_pairing_status", "") for r in candidates if r.get("pdb_pairing_status"))
    binder_type_status_counts = Counter(r.get("binder_type_status", "") for r in candidates if r.get("binder_type_status"))
    review_queue = binder_type_review_queue(candidates)
    qa_counts = Counter(
        flag
        for r in candidates
        for flag in clean(r.get("qa_flags")).split(";")
        if flag
    )
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
- Rows with normalized KD (nM): {n_normalized_kd}/{len(candidates)}
- Rows with antigen UniProt IDs: {n_uniprot}/{len(candidates)}
- Rows with antigen species: {n_antigen_species}/{len(candidates)}
- Rows with binder species hints: {n_binder_species}/{len(candidates)}
- Binder type review queue: {len(review_queue)} rows
- Interaction label: {', '.join(f'{k} {v}' for k, v in label_counts.most_common())}
- Binder types: {', '.join(f'{k} {v}' for k, v in type_counts.most_common())}
- Confidence: {', '.join(f'{k} {v}' for k, v in conf_counts.most_common())}
- Unique antigens: {len(antigen_counts)}
- Binder type status: {', '.join(f'{k} {v}' for k, v in binder_type_status_counts.most_common())}
- PDB pairing status: {', '.join(f'{k} {v}' for k, v in pairing_counts.most_common())}
- Top QA flags: {', '.join(f'{k} {v}' for k, v in qa_counts.most_common(8)) or 'none'}

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
4. Use `scripts/build_patent_review_artifacts.py` to maintain the USPTO /
   Google Patents validation queue for SEQ ID NO, chain role, and claim/example
   context review.
5. Confirm the definition of "designed protein" before splitting ProteinBase /
   de novo binder rows into a separate output module.
6. Promote reviewed rows into `antibody_antigen_curated.csv` and future
   `candidates_v3.csv`.
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
    annotate_pdb_pairing(all_rows)
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

    review_queue = binder_type_review_queue(all_rows)
    write_csv(OUTPUT_DIR / "binder_type_review_queue.csv", REVIEW_QUEUE_COLUMNS, review_queue)
    print(f"  binder_type_review_queue.csv: {len(review_queue)} rows")

    write_memo(all_rows, sample)
    print("\nDone.")


if __name__ == "__main__":
    main()
