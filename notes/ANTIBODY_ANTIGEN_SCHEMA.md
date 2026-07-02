# Antibody-Antigen Schema

## Candidate Table

`outputs/antibody_antigen_candidates.csv` contains sequence-positive candidates
that still need manual review.

Required core columns:

- `binder_name`
- `binder_type`
- `binder_sequence`
- `antigen_name`
- `target_name`
- `target_sequence`
- `interaction_label`
- `affinity_value`
- `affinity_unit`
- `evidence_text`
- `source_reference`
- `confidence`
- `needs_review`

Normalized affinity (added during curation):

- `normalized_KD_nM` — KD converted to nanomolar. Populated during curation after
  verifying that the raw `affinity_value` + `affinity_unit` are consistent with the
  source. `affinity_unit` in the seed is `source_KD` (raw molar values from
  ProteinBase); the curation step should convert to nM and record the result here.

Allowed manual antibody type classes:

- `nanobody`
- `scFv`
- `heavy chain`
- `light chain`
- `CDR-only`
- `other`

Note on `antigen_name` vs `target_name`: in the seed, these are identical for all
rows. They exist as separate fields because during curation the reviewer may
determine that the immunological antigen differs from the interaction target (e.g.,
the binder targets a receptor but was raised against a soluble ectodomain). Update
`antigen_name` when the distinction matters; leave them equal when they are the same
entity.

## Curated Table

`outputs/antibody_antigen_curated.csv` starts as a header-only table. Promote a
row from candidates only after manual review confirms:

- The target is truly the antigen.
- The binder-target label is supported by source evidence.
- The antibody fragment type is correctly classified.
- The sequence-to-antigen pairing is correct.

## Literature Reference Table

`pubmed_references.csv` is an evidence table, not a candidate
sequence-pair table. PubMed rows should be used to support manual curation,
source checking, and later KD/assay extraction. Required columns:

- `antigen_name`
- `query`
- `pmid`
- `title`
- `journal`
- `pubdate`
- `doi`
- `authors`
- `source_reference`
- `evidence_text`
- `url`

## Patent-Scale Additions

Future patent extraction should add source-level sequence listing fields before
curation, including patent identifier, claim/example/table context, `SEQ ID NO`,
heavy/light-chain pairing evidence, and extracted evidence sentence offsets.
