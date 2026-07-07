#!/usr/bin/env python3
"""
Step 2: Normalize antigen names across sources.

Applies a curated mapping to unify variant antigen names (e.g.,
'B-lymphocyte antigen CD20' → 'CD20', 'human tnfa' → 'TNF-alpha').

Input:  any CSV with an 'antigen_name' column
Output: same CSV with antigen_name normalized in-place
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

ANTIGEN_MAP: dict[str, str] = {
    "egfr": "EGFR",
    "b-lymphocyte antigen cd20": "CD20",
    "cd20 antigen": "CD20",
    "ms4a1": "CD20",
    "human tnfa": "TNF-alpha",
    "tnf alpha": "TNF-alpha",
    "tumor necrosis factor": "TNF-alpha",
    "tumor necrosis factor alpha": "TNF-alpha",
    "pd l1": "PD-L1",
    "programmed cell death 1 ligand 1": "PD-L1",
    "programmed death-ligand 1": "PD-L1",
    "cd274": "PD-L1",
    "programmed cell death protein 1": "PD-1",
    "pd 1": "PD-1",
    "pdcd1": "PD-1",
    "epidermal growth factor receptor": "EGFR",
    "erbb1": "EGFR",
    "her1": "EGFR",
    "human epidermal growth factor receptor 2": "HER2",
    "erbb2": "HER2",
    "neu": "HER2",
    "vascular endothelial growth factor a": "VEGF-A",
    "vascular endothelial growth factor": "VEGF-A",
    "vegfa": "VEGF-A",
    "spike glycoprotein": "SARS-CoV-2 Spike",
    "spike protein s1": "SARS-CoV-2 Spike",
    "sars-cov-2 spike protein": "SARS-CoV-2 Spike",
    "spike protein": "SARS-CoV-2 Spike",
    "interleukin-6": "IL-6",
    "il6": "IL-6",
    "il7r": "IL-7R",
    "interleukin-17a": "IL-17A",
    "il17a": "IL-17A",
    "cytotoxic t-lymphocyte protein 4": "CTLA-4",
    "cd152": "CTLA-4",
    "nerve growth factor": "NGF",
    "beta-ngf": "NGF",
    "respiratory syncytial virus fusion glycoprotein": "RSV F protein",
    "rsv f protein": "RSV F protein",
    "fusion glycoprotein f0": "RSV F protein",
    "nipah glycoprotein g": "Nipah glycoprotein G",
    "human serum albumin": "Human serum albumin",
    "rbx1": "RBX1",
    "erbb3": "HER3",
    "human epidermal growth factor receptor 3": "HER3",
    "b-cell maturation antigen": "BCMA",
    "tnfrsf17": "BCMA",
    "prostate-specific membrane antigen": "PSMA",
    "folh1": "PSMA",
    "glutamate carboxypeptidase ii": "PSMA",
    "fibroblast activation protein alpha": "FAP",
    "mesothelin": "Mesothelin",
    "glypican-3": "GPC3",
    "glypican 3": "GPC3",
    "delta-like protein 3": "DLL3",
    "trophoblast cell-surface antigen 2": "TROP2",
    "tacstd2": "TROP2",
    "nectin-4": "Nectin-4",
    "pvrl4": "Nectin-4",
    "claudin-18": "CLDN18.2",
    "claudin 18": "CLDN18.2",
    "epithelial cell adhesion molecule": "EpCAM",
    "cd326": "EpCAM",
    "interleukin-13": "IL-13",
    "il13": "IL-13",
    "interleukin-4 receptor": "IL-4Ra",
    "il4r": "IL-4Ra",
    "interleukin-5": "IL-5",
    "il5": "IL-5",
    "interleukin-23": "IL-23",
    "il23": "IL-23",
    "interleukin-33": "IL-33",
    "il33": "IL-33",
    "interleukin-31": "IL-31",
    "il31": "IL-31",
    "proprotein convertase subtilisin/kexin type 9": "PCSK9",
    "complement c5": "C5",
    "complement component 5": "C5",
    "calcitonin gene-related peptide": "CGRP",
    "angiopoietin-2": "Ang2",
    "angiopoietin 2": "Ang2",
    "angpt2": "Ang2",
    "b-cell activating factor": "BAFF",
    "tnfsf13b": "BAFF",
    "blys": "BAFF",
    "receptor activator of nuclear factor kappa-b ligand": "RANKL",
    "tnfsf11": "RANKL",
    "cd47 antigen": "CD47",
    "integrin-associated protein": "CD47",
    "signal regulatory protein alpha": "SIRPa",
    "sirpa": "SIRPa",
    "cd172a": "SIRPa",
    "t-cell immunoreceptor with ig and itim domains": "TIGIT",
    "hepatitis a virus cellular receptor 2": "TIM-3",
    "havcr2": "TIM-3",
    "lymphocyte activation gene 3": "LAG-3",
    "lag3": "LAG-3",
    "cd223": "LAG-3",
    "b7-h3": "B7-H3",
    "cd276": "B7-H3",
    "thymic stromal lymphopoietin": "TSLP",
    "immunoglobulin e": "IgE",
    "ige": "IgE",
    "kinase insert domain receptor": "VEGFR2",
    "kdr": "VEGFR2",
    "vegfr-2": "VEGFR2",
}


def normalize_antigen(name: str) -> str:
    if not name:
        return name
    key = re.sub(r"\s+", " ", name.strip().lower())
    return ANTIGEN_MAP.get(key, name)


def normalize_file(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "antigen_name" not in (reader.fieldnames or []):
            print(f"  {path.name}: no 'antigen_name' column, skipping")
            return
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    changes = 0
    for row in rows:
        old = row["antigen_name"]
        new = normalize_antigen(old)
        if new != old:
            row["antigen_name"] = new
            changes += 1

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"  {path.name}: {changes} antigen names normalized ({len(rows)} rows)")


def main() -> None:
    if len(sys.argv) < 2:
        project = Path(__file__).resolve().parents[1]
        targets = [
            project / "outputs" / "antibody_antigen_candidates.csv",
            project / "data" / "pdb" / "pdb_ab_ag_pairs.csv",
            project / "data" / "patent" / "patent_ab_ag_pairs.csv",
        ]
        targets = [t for t in targets if t.exists()]
    else:
        targets = [Path(p) for p in sys.argv[1:]]

    if not targets:
        print("No files to normalize.")
        return

    for path in targets:
        normalize_file(path)

    if len(targets) == 1:
        with targets[0].open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        counts = Counter(r.get("antigen_name", "") for r in rows)
        top = counts.most_common(10)
        print(f"\nTop 10 antigens: {', '.join(f'{n}({c})' for n, c in top)}")


if __name__ == "__main__":
    main()
