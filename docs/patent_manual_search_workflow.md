# Patent Manual Search Workflow

This note is written for live demonstration. The goal is to show that the
pipeline is not a black box: the same patent sequence evidence can be searched
manually in NCBI, USPTO, and Google Patents, then scaled with code.

## Why three patent interfaces are used

| Interface | Best use in this project | What it gives us |
|-----------|--------------------------|------------------|
| NCBI Protein, PAT division | Fast retrieval of patent-deposited protein FASTA records | Accession, protein title, amino-acid sequence |
| USPTO Patent Public Search | Official US patent / application context | Patent title, claims, description, sequence-listing references |
| USPTO sequence listing site | Sequence listings for granted / published US patents | Tables, sequence listings, mega sequence items |
| Google Patents | Fast manual cross-check across patent families | Patent family, assignee, claims, full-text search, "SEQ ID NO" context |

The automated dataset currently uses NCBI Protein PAT for sequence retrieval
because it exposes FASTA records through E-utilities. USPTO and Google Patents
are added as manual validation routes because they are better for checking
whether a sequence belongs to a patent claim, example, heavy chain, light chain,
or antigen-binding region.

## Demo 1: NCBI Protein PAT

**Website:** https://www.ncbi.nlm.nih.gov/protein/

Use this when the professor asks: "NCBI 你到哪个网站怎么搜?"

### Example target: PD-1

Paste this query into the NCBI Protein search box:

```text
(anti-PD-1 OR anti-PDCD1 OR nivolumab OR pembrolizumab OR cemiplimab)
AND (antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody OR Fab OR heavy chain OR light chain OR CDR OR variable region)
AND "PAT"[Division]
```

What to show:

1. The search is in **NCBI Protein**, not PubMed.
2. `"PAT"[Division]` restricts the results to GenBank patent protein records.
3. Open one result and show that the record has a protein accession and amino-acid sequence.
4. Click **FASTA** to show that the sequence is directly downloadable.
5. Explain that `scripts/fetch_patent_sequences.py` performs the same search by API:
   `esearch` gets record IDs, then `efetch` downloads FASTA.

### Example target: EGFR

```text
(anti-EGFR OR anti-ERBB1 OR cetuximab OR panitumumab OR necitumumab)
AND (antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody OR Fab OR heavy chain OR light chain OR CDR OR variable region)
AND "PAT"[Division]
```

### Example target: BCMA

```text
(anti-BCMA OR anti-TNFRSF17 OR teclistamab OR elranatamab)
AND (antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody OR Fab OR heavy chain OR light chain OR CDR OR variable region)
AND "PAT"[Division]
```

## Demo 2: USPTO Patent Public Search

**Website:** https://www.uspto.gov/patents/search/patent-public-search

USPTO Patent Public Search is the official web app for US patents and patent
application publications. Use this to show patent context after NCBI has found
a sequence or after Google Patents has found a family.

### Basic search example

Open Patent Public Search, choose **Basic Search**, then search:

```text
"PD-1" antibody variable region
```

or:

```text
pembrolizumab antibody sequence
```

What to show:

1. Results are patents / published applications, not protein FASTA records.
2. Open a result and search within the document for `SEQ ID NO`.
3. Show that patent documents often describe sequence IDs, heavy chain,
   light chain, CDRs, claims, and examples.
4. Explain why this is the next curation step: NCBI gives sequences, but USPTO
   helps validate which sequence is heavy chain, light chain, CDR, or claimed
   antibody.

## Demo 3: USPTO Issued and Published Sequences

**Entry page:** https://www.uspto.gov/patents/search

The USPTO search page links to **Issued and Published Sequences**, which is
used when we already know the patent number or publication number and want the
official sequence listing.

What to show:

1. Start from USPTO patent search page.
2. Open **Issued and Published Sequences**.
3. Search by a known patent or publication number found from Patent Public
   Search / Google Patents.
4. Use this as the source for future SEQ ID NO parsing.

## Demo 4: Google Patents

**Website:** https://patents.google.com/

Google Patents is useful for fast manual checking because it supports phrase
search, metadata restrictions, fielded search, and patent-family navigation.

### Example query

```text
("PD-1" OR PDCD1 OR pembrolizumab OR nivolumab) (antibody OR immunoglobulin OR scFv OR Fab OR nanobody) country:US
```

### Example fielded query

```text
TI=("PD-1" antibody) OR AB=("PD-1" antibody) OR CL=("PD-1" antibody)
```

What to show:

1. Open a patent result.
2. Search within the page for `SEQ ID NO`.
3. Show claims / description / sequence listing references.
4. Note the assignee and patent family.
5. If a sequence record from NCBI has a patent accession or title, use Google
   Patents to find the corresponding patent family and context.

## How manual search connects to the code

```text
Manual query design
      |
      v
Target keyword table
      |
      v
NCBI ESearch: find PAT protein record IDs
      |
      v
NCBI EFetch: download FASTA sequences
      |
      v
Length / amino-acid filters
      |
      v
UniProt antigen sequence attachment
      |
      v
Merge, deduplicate, QA flags
      |
      v
Manual USPTO / Google Patents validation for SEQ ID NO and chain pairing
```

## Current limitation to state honestly

The current `NCBI Patent` rows are **target-annotated sequence candidates**.
They are not yet fully curated patent antibody pairs. The next step is to parse
patent full text / sequence listings to recover:

- `SEQ ID NO`
- heavy-chain / light-chain pairing
- CDR boundaries
- whether the sequence is claimed, example-only, or background prior art
- the exact antigen-binding context

That is why `outputs/manual_validation_sample.csv` and
`outputs/binder_type_review_queue.csv` exist.

