# Project Explanation for Presenter

This file is written in plain language so the project can be explained during a
meeting without getting lost in code details.

## What the project is trying to build

The project is building a database of antibody / binder sequences paired with
their antigen / target sequences.

One row means:

```text
binder sequence -> antigen sequence | source | evidence | QA status
```

This is different from a PubMed-only literature database. PubMed can tell us
which papers discuss an antibody, but the database needs amino-acid sequences.
That is why the main sequence sources are PDB and NCBI Protein PAT.

## Current implemented data sources

| Source | What it contributes | Current role |
|--------|---------------------|--------------|
| ProteinBase seed | 300 seed antibody-like binder-antigen rows | v1 seed and comparison baseline |
| RCSB PDB | co-crystal structure-derived antibody-antigen protein chains | high-confidence structural candidates |
| NCBI Protein PAT | patent-deposited antibody-like protein sequences | large-scale patent sequence expansion |
| UniProt | canonical antigen sequences for searched targets | attaches target sequences to patent rows |
| PubMed | literature metadata by antigen | supporting evidence layer, not direct pair rows |

## Current scale

- 3,733 total sequence-pair candidates.
- 2,679 final NCBI Patent rows after merge/deduplication.
- 754 PDB rows.
- 300 ProteinBase seed rows.
- 325 unique antigen names in the current master table.
- 120 PubMed reference rows for literature support.

## How the code works

### Step 1: Choose target keywords

The project defines therapeutic targets such as PD-1, EGFR, BCMA, TROP2,
IL-17A, PCSK9, RSV F, and others. Each target has synonyms:

```text
PD-1 = anti-PD-1 OR anti-PDCD1 OR nivolumab OR pembrolizumab OR cemiplimab
```

The search also adds antibody-related terms:

```text
antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody OR Fab
OR heavy chain OR light chain OR CDR OR variable region
```

### Step 2: Search NCBI Protein PAT

The code searches NCBI Protein, restricted to patent records:

```text
({target synonyms})
AND ({antibody terms})
AND "PAT"[Division]
```

In code this is done by `scripts/fetch_patent_sequences.py`:

1. `esearch` finds NCBI protein record IDs.
2. `efetch` downloads FASTA sequences.
3. The script filters sequences by length and amino-acid characters.
4. Each patent sequence is paired with the target's UniProt antigen sequence.

### Step 3: Search PDB structures

`scripts/fetch_pdb_complexes.py` searches RCSB PDB for antibody-antigen
co-crystal entries. It downloads polymer entity metadata and sequences, then
pairs antibody-like chains with antigen-like chains inside the same PDB complex.

### Step 4: Normalize and merge

`scripts/merge_candidates.py` combines all sources into one schema.

Important cleanup:

- computes `sequence_pair_key` from binder and antigen sequence;
- removes duplicate sequence pairs;
- normalizes antigen names;
- fills KD in nM when possible;
- flags ambiguous residues;
- marks rows needing manual review.

### Step 5: Serve the app

`app.py` reads `antibody_antigen_master_v2.csv` directly. The app is not the
data generator. It is the browsing and review layer.

## Where LLM was used

The LLM was used for workflow design, keyword expansion, documentation, and QA
planning. It was **not** used to invent sequence data.

Example prompt used for keyword expansion:

```text
For each therapeutic antibody target, propose search synonyms for NCBI Protein
PAT retrieval. Include gene symbols, CD names, target aliases, and approved or
late-stage antibody drug names. Return a conservative boolean query. Avoid
generic terms that would retrieve unrelated proteins.
```

Example prompt used for QA planning:

```text
Given the current antibody-antigen master table columns and source types,
identify the highest-risk rows for manual validation. Prioritize patent rows
where binder type is ambiguous, sequence length is unusual, or the row lacks
chain-pairing evidence.
```

The important rule is: the LLM can propose search terms and review priorities,
but every database row must still come from a traceable source such as NCBI,
PDB, UniProt, ProteinBase, PubMed, USPTO, or Google Patents.

## What was improved after the meeting

| Problem raised in meeting | Change made |
|---------------------------|-------------|
| Data volume was too small | NCBI target list expanded to 64 targets and retmax raised to 200 |
| Zora could not manually demonstrate search | Added NCBI / USPTO / Google Patents manual workflow |
| NCBI search was not broad enough | Added target synonym table with hit/kept counts |
| Natural and designed protein were mixed conceptually | Slides now show two-track architecture |
| Need real framework and real evidence | Slides use actual master-table counts, charts, code snippets, and search examples |

## Main limitations to say honestly

1. NCBI PAT rows are candidates, not fully curated antibody pairs yet.
2. Patent rows still need SEQ ID NO extraction and heavy/light chain pairing.
3. PDB rows need additional literature mining to fill KD values.
4. Designed protein needs a stricter definition before splitting ProteinBase rows.
5. The next high-value version is `v3`: manually curated, patent-context-backed,
   chain-paired entries.
