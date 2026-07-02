# Antibody-Antigen PairingDB

**A searchable antibody-antigen sequence-pair candidate database, built as the
protein-interaction extension of CellReprogramDB.**

This project follows the same architecture as
`/Users/chen/reprogramming_project`: source fetch scripts and curation outputs
produce one root master table, and the Streamlit app reads that master table
directly.

## Overview

Antibody-antigen pairing requires two related but different evidence layers:

```text
antibody / binder sequence -> antigen / target sequence | source | evidence
```

The current database contains **1,310 sequence-pair candidates** from three
structured sources and **120 PubMed reference rows** as a literature evidence
layer. PubMed rows are not counted as sequence-pair candidates unless a paper or
supplement is parsed to recover paired amino-acid sequences.

## Pipeline

```text
ProteinBase seed / RCSB PDB / NCBI GenBank PAT
      |
      v
scripts/fetch_pdb_complexes.py       -> data/pdb/pdb_ab_ag_pairs.csv
scripts/fetch_patent_sequences.py    -> data/patent/patent_ab_ag_pairs.csv
scripts/fetch_pubmed_references.py   -> pubmed_references.csv
      |
      v
scripts/normalize_antigens.py        -> normalized antigen names
      |
      v
scripts/merge_candidates.py          -> antibody_antigen_master_v2.csv
                                      -> outputs/candidates_v1.csv
                                      -> outputs/candidates_v2.csv
                                      -> outputs/manual_validation_sample.csv
                                      -> outputs/validation_memo.md
      |
      v
app.py                               -> Streamlit web app
```

## Data Summary

- **1,310** antibody-antigen sequence-pair candidates
- **283** unique antigens
- **3** sequence sources: RCSB PDB (754), ProteinBase seed (300), NCBI Patent (256)
- **7** binder types: nanobody, heavy chain, light chain, scFv, Fab, other, antibody fragment
- **90** rows with affinity values
- **120** PubMed references in `pubmed_references.csv`

## Versions

- `outputs/candidates_v1.csv`: 300 rows, ProteinBase seed only
- `outputs/candidates_v2.csv`: 1,310 rows, ProteinBase + PDB + NCBI Patent
- `antibody_antigen_master_v2.csv`: current app-facing master table
- `outputs/antibody_antigen_curated.csv`: header-only curated table until manual review
- `candidates_v3.csv`: planned manually curated, literature-supported subset

## Repository Structure

```text
app.py                              Streamlit web application
antibody_antigen_master_v2.csv      App-facing master candidate table
pubmed_references.csv               PubMed metadata for curation evidence
requirements.txt                    Streamlit app dependencies
requirements-pipeline.txt           Fetch / pipeline dependencies

scripts/
  fetch_pdb_complexes.py            RCSB PDB complex retrieval
  fetch_patent_sequences.py         NCBI GenBank PAT + UniProt retrieval
  fetch_pubmed_references.py        PubMed literature reference retrieval
  normalize_antigens.py             Antigen-name normalization
  merge_candidates.py               Versioned merge and master-table writer
  check_project.py                  Local consistency check
  build_antibody_antigen_candidates.py  Legacy seed-only builder

data/
  source_pilots/                    Protein interaction pilot seed files
  pdb/                              PDB intermediate CSV
  patent/                           Patent intermediate CSV
  literature/                       Historical PubMed reference copy
  raw/                              Raw pilot source CSVs

outputs/
  candidates_v1.csv
  candidates_v2.csv
  antibody_antigen_candidates.csv
  antibody_antigen_curated.csv
  manual_validation_sample.csv
  validation_memo.md

notes/
  ANTIBODY_ANTIGEN_SCHEMA.md
  PROTEIN_INTERACTION_PIPELINE_source.md
```

## Setup

Install only the app dependencies to browse the database:

```bash
pip install -r requirements.txt
streamlit run app.py
```

Install pipeline dependencies only when refreshing source data:

```bash
pip install -r requirements-pipeline.txt
```

## Run The Pipeline

The checked-in master table can be browsed without refreshing remote sources.
To rebuild locally:

```bash
python3 scripts/fetch_pdb_complexes.py
python3 scripts/fetch_patent_sequences.py
python3 scripts/fetch_pubmed_references.py --max-targets 30 --retmax 5
python3 scripts/normalize_antigens.py
python3 scripts/merge_candidates.py
python3 scripts/check_project.py
```

## Streamlit Cloud Deployment

This project should deploy the same way as CellReprogramDB:

1. Commit `app.py`, `requirements.txt`, `antibody_antigen_master_v2.csv`,
   `pubmed_references.csv`, source scripts, and documentation to GitHub.
2. Create a Streamlit Community Cloud app from the GitHub repository.
3. Set the entry point to `app.py`.
4. Let Streamlit install `requirements.txt` and run the app.

Pipeline files and larger source/intermediate CSVs can stay in the repo while
the app continues to read the root master table.

## Quality Notes

- PDB-derived pairs are high-confidence structure candidates but still need
  manual review. A small number of possible antibody-chain-as-antigen cases are
  flagged in `qa_flags`.
- NCBI Patent rows are target-annotated sequence candidates. They need patent
  context review for heavy/light pairing, `SEQ ID NO`, and claim/example support.
- PubMed references support curation but are not direct sequence-pair evidence.

See `notes/ANTIBODY_ANTIGEN_SCHEMA.md` for the schema.
