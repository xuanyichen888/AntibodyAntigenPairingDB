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
- `antigen_uniprot_id`
- `binder_species`
- `antigen_species`
- `interaction_label`
- `affinity_value`
- `affinity_unit`
- `normalized_KD_nM`
- `evidence_text`
- `source_reference`
- `confidence`
- `needs_review`
- `complex_group_id`
- `pdb_pairing_status`
- `paired_chain_roles`
- `binder_type_status`
- `binder_type_suggestion`
- `qa_flags`

Normalized affinity:

- `normalized_KD_nM` — KD converted to nanomolar. The merge step converts
  ProteinBase `source_KD` molar values into this field. Curators should verify
  the source assay before promoting rows to the curated table.

Source and QA fields:

- `antigen_uniprot_id` — UniProt accession for the antigen when available.
  Patent rows are backfilled from the target map used by
  `scripts/fetch_patent_sequences.py`.
- `qa_flags` — semicolon-separated machine-readable warnings. Current merge-time
  flags include `binder_sequence_contains_x`, `target_sequence_contains_x`,
  `binder_sequence_invalid_chars`, `target_sequence_invalid_chars`,
  `binder_sequence_ambiguous_aa`, `target_sequence_ambiguous_aa`, and
  `possible_pdb_antibody_chain_as_antigen`.
- `binder_species` / `antigen_species` — conservative organism hints inferred
  from UniProt target maps, antigen names, and explicit source text. Blank means
  the pipeline did not have enough evidence.
- `complex_group_id`, `paired_chain_roles`, `pdb_pairing_status` — PDB-only
  complex-level pairing hints. Because the current PDB intermediate table does
  not retain chain IDs, these fields group by PDB ID and antigen rather than
  claiming exact heavy/light chain identity.
- `binder_type_status` — `classified`, `needs_manual_review`, or
  `non_antibody_candidate`. Rows that need review are written to
  `outputs/binder_type_review_queue.csv`.

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
