# Revised Presentation Outline

This outline is the corrected version for the next deck / HTML report. It is
designed to be defensible under professor questioning: every numeric claim is
backed by a checked file, and every unfinished item is explicitly labeled as a
candidate or validation queue.

## Main Framing

**One-sentence message:**  
I reused the CellReprogramDB pipeline pattern, but changed the biological object
from a text-mined reprogramming recipe to a sequence-pair candidate:

```text
binder / antibody sequence -> antigen / target sequence | source | evidence | QA status
```

**Important wording:**  
This is a **v2 candidate database**, not a curated final database. All rows are
currently `candidate_unreviewed`, and v3 will be the curated patent/context
validated subset.

## Slide 1: What changed after the meeting

Purpose: show that the feedback was converted into build tasks.

Real evidence:

- Total candidate rows: 1,310 -> 3,733.
- NCBI Patent final rows: 256 -> 2,679.
- NCBI target search panel: 64 targets.
- New validation outputs:
  - `outputs/patent_search_audit.csv`
  - `outputs/patent_validation_queue.csv`
  - `outputs/patent_validation_patent_summary.csv`

## Slide 2: What this database stores

Purpose: make the project concept simple before code.

Show:

```text
binder sequence -> antigen sequence | source | confidence | QA flags
```

Contrast:

| CellReprogramDB | Antibody-Antigen PairingDB |
|-----------------|----------------------------|
| source cell -> target cell + factors | binder sequence -> antigen sequence |
| PubMed text extraction | PDB / NCBI PAT / ProteinBase sequence evidence |
| recipes_master_v2.csv | antibody_antigen_master_v2.csv |

## Slide 3: Reused CellReprogramDB architecture

Purpose: answer "why this architecture".

Flow:

```text
source-specific fetchers
  -> normalized schema
  -> merge + dedupe + QA
  -> versioned master table
  -> Streamlit app
  -> manual curation queue
```

Say explicitly: I did not create a new standard; I copied the stable
CellReprogramDB pattern and adapted it to sequence-pair evidence.

## Slide 4: Current v2 data status

Purpose: give honest numbers.

Use:

- `antibody_antigen_master_v2.csv`: 3,733 rows.
- `source_type`:
  - NCBI Patent: 2,679
  - PDB: 754
  - ProteinBase: 300
- Unique antigen names: 325.
- PubMed references: 120, evidence layer only.

Do not say "complete" or "curated". Say "expanded candidate set".

## Slide 5: Pipeline with real files

Purpose: show exactly how code runs.

Flow:

```text
scripts/fetch_pdb_complexes.py       -> data/pdb/pdb_ab_ag_pairs.csv
scripts/fetch_patent_sequences.py    -> data/patent/patent_ab_ag_pairs.csv
scripts/fetch_pubmed_references.py   -> pubmed_references.csv
scripts/merge_candidates.py          -> antibody_antigen_master_v2.csv
scripts/build_patent_review_artifacts.py
                                     -> patent validation outputs
app.py                               -> Streamlit browser
```

Use a small code screenshot from `fetch_patent_sequences.py` and
`build_patent_review_artifacts.py`.

## Slide 6: NCBI Protein PAT search logic

Purpose: answer "NCBI 怎么搜".

Example query:

```text
(anti-PD-1 OR anti-PDCD1 OR nivolumab OR pembrolizumab OR cemiplimab)
AND (antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody OR Fab OR heavy chain OR light chain OR CDR OR variable region)
AND "PAT"[Division]
```

Manual demo:

1. Open NCBI Protein.
2. Paste query.
3. Open a result.
4. Click FASTA.
5. Show code repeats the same search for 64 targets.

## Slide 7: Keyword search strategy

Purpose: prove search was widened systematically.

Show 6-8 representative rows from `outputs/patent_search_audit.csv`:

- target
- category
- target query
- UniProt ID
- NCBI hits
- final master rows
- capped at retmax or not

Important explanation:

When `ncbi_hits_reported = 200`, the result may be capped by `retmax=200`.
That means it is a coverage flag, not proof the target is exhausted.

## Slide 8: Results and deduplication

Purpose: show improvement with honest mechanics.

Numbers:

- Raw sources before dedupe: 4,883.
- Removed duplicate sequence pairs: 1,150.
- Final master rows: 3,733.
- Patent intermediate: 3,829.
- Patent final after merge/dedupe: 2,679.

Explain:

Deduplication uses MD5 of binder sequence + target sequence, so repeated
accessions or repeated target-pair rows collapse into one candidate row.

## Slide 9: QA reality

Purpose: avoid overclaiming and show curation maturity.

Numbers:

- All rows are `candidate_unreviewed`.
- All rows have `needs_review=True`.
- Binder type review queue: 2,205 rows.
- Patent validation queue: 2,679 rows.
- Patent-level summary: 172 patents.

Message:

The project is now good enough to search and review systematically. It is not
yet a manually curated antibody-pair gold standard.

## Slide 10: USPTO / Google Patents validation workflow

Purpose: answer "专利网站怎么搜".

Show:

- `outputs/patent_validation_queue.csv`
- `outputs/patent_validation_patent_summary.csv`

Validation fields:

- patent number
- NCBI accession
- inferred Sequence label
- Google Patents URL
- USPTO search URL
- manual patent family verified
- manual SEQ ID NO verified
- manual chain role
- manual claim context

Key wording:

NCBI gives sequence candidates. USPTO / Google Patents validate what those
sequences mean inside the patent document.

## Slide 11: Streamlit app as review interface

Purpose: show this is usable.

Current app tabs:

- Data Table
- Charts
- Pair Detail
- Literature
- Validation

The new Validation tab exposes search audit, patent queue, and patent summary
without requiring the professor to inspect CSV files manually.

## Slide 12: Next decisions

Purpose: ask for professor guidance.

Decision points:

1. What exactly counts as "designed protein"?
2. Should ProteinBase scFv/nanobody rows be natural antibody, designed protein,
   or mixed evidence?
3. For v3, should priority be:
   - sequence-pair confidence,
   - patent claim support,
   - KD extraction,
   - or broader coverage?

Recommendation:

Keep one repository for now, split outputs into two modules after the definition
is confirmed:

```text
natural_antibody_master.csv
designed_protein_master.csv
```

