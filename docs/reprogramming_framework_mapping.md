# Mapping the CellReprogramDB Framework to Antibody-Antigen PairingDB

The antigen-antibody project is intentionally built as a reusable extension of
`/Users/chen/reprogramming_project`. The stable idea is the same: collect
source evidence, normalize it into one schema, version the outputs, run QA, and
serve a Streamlit app from the master table.

## Shared architecture

| CellReprogramDB | Antibody-Antigen PairingDB | Reused idea |
|-----------------|----------------------------|-------------|
| PubMed keyword search | PDB / NCBI PAT / PubMed / ProteinBase source search | source-specific fetchers |
| Abstract / full-text extraction | FASTA / polymer entity / patent sequence extraction | structured evidence extraction |
| recipe schema | sequence-pair schema | one normalized row format |
| source cell, target cell, factors | binder sequence, antigen sequence, source evidence | biological object normalization |
| duplicate recipe marking | MD5 sequence-pair deduplication | stable duplicate keys |
| confidence / needs review | confidence / needs review / QA flags | curation state tracking |
| recipes_master_v2.csv | antibody_antigen_master_v2.csv | app-facing master table |
| Streamlit app | Streamlit app | searchable browser + charts |
| manual QA samples | manual_validation_sample.csv | human curation queue |
| supplementary keyword tables | docs/search_strategy.md | reproducibility package |

## What was implemented in this project

1. **Source fetchers**
   - `scripts/fetch_pdb_complexes.py`: RCSB PDB structure-derived antibody-antigen complexes.
   - `scripts/fetch_patent_sequences.py`: NCBI Protein PAT antibody-like sequences plus UniProt antigen sequences.
   - `scripts/fetch_pubmed_references.py`: PubMed literature metadata as an evidence layer.

2. **Normalized master schema**
   - One row represents `binder / antibody sequence -> antigen / target sequence`.
   - Common columns include sequence fields, source metadata, confidence,
     affinity, QA flags, duplicate key, and curation status.

3. **Versioning**
   - `outputs/candidates_v1.csv`: ProteinBase seed only.
   - `outputs/candidates_v2.csv`: expanded PDB + NCBI PAT + ProteinBase set.
   - Future `candidates_v3.csv`: manually curated, literature/patent validated subset.

4. **QA and curation**
   - Sequence-pair MD5 deduplication.
   - Ambiguous residue flags.
   - KD normalization into `normalized_KD_nM`.
   - PDB complex-level pairing metadata.
   - Patent binder-type review queue.

5. **App layer**
   - `app.py` reads the master table directly.
   - Tabs: Data Table, Charts, Pair Detail, Literature.
   - Filters: source, binder type, antigen, confidence, binder type status,
     interaction label, sequence length.

## What is different from CellReprogramDB

| Difference | Why it matters |
|------------|----------------|
| The biological unit is a sequence pair, not a text recipe | We need exact amino-acid strings, not only natural-language extraction |
| Patent data has many sequence fragments | Heavy/light pairing and SEQ ID NO validation become a separate curation problem |
| PDB data is structural but incomplete for affinity | It gives high-confidence physical complexes, but KD often needs literature backfill |
| PubMed is evidence support, not direct sequence source | PubMed abstracts rarely contain full antibody sequences |
| Designed proteins need a separate track | De novo/miniprotein binders should not be mixed with natural antibody evidence |

## Two-subproject framing for the next version

The code can stay in one repository while the data model is split into two
tracks:

```text
Antibody-Antigen PairingDB
      |
      +-- Track A: Natural Antibody
      |       Sources: PDB, NCBI GenBank PAT, PubMed/PMC, UniProt
      |       Output: natural_antibody_master.csv
      |
      +-- Track B: Designed Protein
              Sources: ProteinBase seed, de novo binder datasets, RFdiffusion/ProteinMPNN papers
              Output: designed_protein_master.csv
```

This avoids prematurely classifying all ProteinBase rows as designed protein.
For the presentation, the split should be shown as the target architecture. The
current implemented master table is still the shared v2 candidate table.

