# Search Strategy — Antibody-Antigen Sequence Pair Collection

This document records the systematic search strategy used to construct the
Antibody-Antigen PairingDB. It is intended as supplementary material to
demonstrate reproducibility.

## Data Sources

| Source | Database | API / Interface | Records |
|--------|----------|-----------------|---------|
| RCSB PDB | Protein Data Bank | RCSB Search API v2 + REST API | see below |
| NCBI GenBank PAT | NCBI Protein (patent division) | E-utilities (esearch + efetch) | see below |
| ProteinBase | Literature-curated seed | Manual CSV | 300 |

---

## Source 1: RCSB PDB Co-crystal Structures

**Database:** RCSB Protein Data Bank (https://www.rcsb.org)

**Search strategy:** Full-text search for antibody-antigen co-crystal structures
using the RCSB Search API v2. Each hit is expanded via the REST API to extract
all polymer entities, then antibody and antigen chains are paired by PDB complex.

**Query:** `full_text: "antibody antigen"` via RCSB Search API

**Filters applied:**
- Polymer entity type: protein
- Sequence extraction: all polymer entities per PDB entry
- Pairing: antibody and antigen chains grouped by PDB ID

**Post-processing:**
- Binder type classified from entity description (VH, VL, VHH, scFv, Fab, nanobody)
- Complex group ID assigned for VH/VL chain pairing within same PDB entry
- Sequences with non-standard residues (X) flagged in qa_flags

**Script:** `scripts/fetch_pdb_complexes.py`

---

## Source 2: NCBI GenBank PAT (Patent-deposited Sequences)

**Database:** NCBI Protein — PAT division (https://www.ncbi.nlm.nih.gov/protein)

**Search strategy:** For each therapeutic target, a combined query searches for
patent-deposited protein sequences that mention both the target and antibody-related
terms. Antigen sequences are retrieved from UniProt for each target.

**Query template:**
```
({target_synonyms}) AND (antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv
OR nanobody OR Fab OR heavy chain OR light chain OR CDR OR variable region)
AND "PAT"[Division]
```

**Parameters:**
- `db`: protein
- `retmax`: 200 per target
- `retmode`: json (esearch), text/fasta (efetch)

**Sequence filters:**
- Length: 80–600 amino acids (variable domain range)
- Composition: standard amino acid characters only (ACDEFGHIKLMNPQRSTVWYXBZU)
- Deduplication by NCBI accession

**Antigen sequences:** Retrieved from UniProt REST API using curated UniProt IDs

### Target Keyword Table

| Target | Category | Query Synonyms | UniProt ID | NCBI Hits | Kept |
|--------|----------|---------------|------------|-----------|------|
| PD-1 | Checkpoint | anti-PD-1, anti-PDCD1, nivolumab, pembrolizumab, cemiplimab | Q15116 | 200 | 137 |
| PD-L1 | Checkpoint | anti-PD-L1, anti-CD274, atezolizumab, durvalumab, avelumab | Q9NZQ7 | 200 | 63 |
| CTLA-4 | Checkpoint | anti-CTLA-4, anti-CTLA4, ipilimumab, tremelimumab | P16410 | 200 | 137 |
| TIGIT | Checkpoint | anti-TIGIT, tiragolumab, vibostolimab, domvanalimab | Q495A1 | 200 | 179 |
| LAG-3 | Checkpoint | anti-LAG-3, anti-LAG3, relatlimab, favezelimab | P18627 | 200 | 97 |
| TIM-3 | Checkpoint | anti-TIM-3, anti-TIM3, anti-HAVCR2, cobolimab, sabatolimab | Q8TDQ0 | 193 | 103 |
| VISTA | Checkpoint | anti-VISTA, anti-VSIR, anti-B7-H5 | Q9H7M9 | 46 | 18 |
| BTLA | Checkpoint | anti-BTLA, anti-CD272 | Q7Z6A9 | 0 | 0 |
| CD47 | Checkpoint | anti-CD47, magrolimab, lemzoparlimab | Q08722 | 200 | 97 |
| SIRPa | Checkpoint | anti-SIRPa, anti-SIRPA, anti-CD172a | P78324 | 0 | 0 |
| HER2 | Oncology | anti-HER2, anti-ERBB2, trastuzumab, pertuzumab, margetuximab | P04626 | 200 | 129 |
| EGFR | Oncology | anti-EGFR, anti-ERBB1, cetuximab, panitumumab, necitumumab | P00533 | 200 | 98 |
| HER3 | Oncology | anti-HER3, anti-ERBB3, patritumab | P21860 | 200 | 52 |
| CD20 | Oncology | anti-CD20, anti-MS4A1, rituximab, obinutuzumab, ofatumumab | P11836 | 200 | 99 |
| CD19 | Oncology | anti-CD19, blinatumomab, loncastuximab, tafasitamab | P15391 | 200 | 105 |
| CD3 | Oncology | anti-CD3, anti-CD3E, blinatumomab, teplizumab | P07766 | 200 | 60 |
| CD38 | Oncology | anti-CD38, daratumumab, isatuximab, mezagitamab | P28907 | 178 | 103 |
| BCMA | Oncology | anti-BCMA, anti-TNFRSF17, teclistamab, elranatamab | Q02223 | 200 | 119 |
| GPRC5D | Oncology | anti-GPRC5D, talquetamab | Q9NZD1 | 2 | 2 |
| FcRH5 | Oncology | anti-FcRH5, anti-FCRL5, cevostamab | Q96RD9 | 0 | 0 |
| CD33 | Oncology | anti-CD33, gemtuzumab | P20138 | 32 | 25 |
| CD22 | Oncology | anti-CD22, inotuzumab, moxetumomab | P20273 | 200 | 16 |
| CD30 | Oncology | anti-CD30, anti-TNFRSF8, brentuximab | P28908 | 200 | 80 |
| TROP2 | Oncology | anti-TROP2, anti-TACSTD2, sacituzumab, datopotamab | P09758 | 97 | 56 |
| Nectin-4 | Oncology | anti-Nectin-4, anti-PVRL4, enfortumab | Q96NY8 | 10 | 5 |
| CLDN18.2 | Oncology | anti-claudin 18, anti-CLDN18, zolbetuximab | P56856 | 17 | 17 |
| DLL3 | Oncology | anti-DLL3, rovalpituzumab, tarlatamab | Q9NYJ7 | 200 | 133 |
| Mesothelin | Oncology | anti-mesothelin, anti-MSLN, anetumab | Q13421 | 160 | 59 |
| MET | Oncology | anti-MET, anti-HGFR, onartuzumab, amivantamab | P08581 | 200 | 40 |
| PSMA | Oncology | anti-PSMA, anti-FOLH1, capromab | Q04609 | 200 | 188 |
| FAP | Oncology | anti-FAP, anti-fibroblast activation protein | Q12884 | 10 | 6 |
| GPC3 | Oncology | anti-GPC3, anti-glypican-3, codrituzumab | P51654 | 200 | 109 |
| B7-H3 | Oncology | anti-B7-H3, anti-CD276, enoblituzumab, omburtamab | Q5ZPR3 | 200 | 106 |
| EpCAM | Oncology | anti-EpCAM, anti-CD326, catumaxomab, solitomab | P16422 | 76 | 44 |
| TNF-alpha | Cytokine | anti-TNF, anti-TNFA, adalimumab, infliximab, golimumab, certolizumab | P01375 | 200 | 61 |
| IL-6 | Cytokine | anti-IL-6, anti-IL6, siltuximab, olokizumab | P05231 | 200 | 77 |
| IL-6R | Cytokine | anti-IL-6R, anti-IL6R, tocilizumab, sarilumab | P08887 | 200 | 75 |
| IL-17A | Cytokine | anti-IL-17A, anti-IL17A, secukinumab, ixekizumab, bimekizumab | Q16552 | 77 | 34 |
| IL-17RA | Cytokine | anti-IL-17RA, anti-IL17RA, brodalumab | Q96F46 | 7 | 2 |
| IL-4Ra | Cytokine | anti-IL-4R, anti-IL4R, dupilumab | P24394 | 76 | 43 |
| IL-13 | Cytokine | anti-IL-13, anti-IL13, tralokinumab, cendakimab | P35225 | 200 | 65 |
| IL-5 | Cytokine | anti-IL-5, anti-IL5, mepolizumab, reslizumab | P05113 | 30 | 13 |
| IL-5Ra | Cytokine | anti-IL-5R, anti-IL5RA, benralizumab | Q01344 | 7 | 2 |
| IL-23 | Cytokine | anti-IL-23, anti-IL23, guselkumab, risankizumab, tildrakizumab | Q9NPF7 | 200 | 101 |
| IL-12/23 p40 | Cytokine | anti-IL-12, anti-p40, ustekinumab | P29460 | 126 | 43 |
| IL-33 | Cytokine | anti-IL-33, anti-IL33, itepekimab, tozorakimab | O95760 | 200 | 113 |
| IL-31 | Cytokine | anti-IL-31, anti-IL31, nemolizumab | Q6EBC2 | 200 | 88 |
| TSLP | Cytokine | anti-TSLP, tezepelumab | Q969D9 | 163 | 65 |
| IgE | Cytokine | anti-IgE, omalizumab, ligelizumab | P01854 | 43 | 20 |
| VEGF | Soluble | anti-VEGF, anti-VEGFA, bevacizumab, ranibizumab, brolucizumab | P15692 | 200 | 116 |
| VEGFR2 | Soluble | anti-VEGFR2, anti-KDR, ramucirumab | P35968 | 9 | 4 |
| NGF | Soluble | anti-NGF, tanezumab, fasinumab, fulranumab | P01138 | 200 | 85 |
| CGRP | Soluble | anti-CGRP, fremanezumab, galcanezumab, eptinezumab | P06881 | 200 | 62 |
| CGRPr | Soluble | anti-CGRP receptor, erenumab, anti-CALCRL | Q16602 | 7 | 2 |
| PCSK9 | Soluble | anti-PCSK9, alirocumab, evolocumab | Q8NBP7 | 200 | 61 |
| C5 | Soluble | anti-C5, eculizumab, ravulizumab | P01031 | 183 | 77 |
| RANKL | Soluble | anti-RANKL, anti-TNFSF11, denosumab | O14788 | 7 | 2 |
| Ang2 | Soluble | anti-angiopoietin-2, anti-ANGPT2, faricimab | O15123 | 51 | 46 |
| BAFF | Soluble | anti-BAFF, anti-TNFSF13B, belimumab | Q9Y275 | 7 | 2 |
| RSV F | Infectious | anti-RSV, anti-respiratory syncytial virus, palivizumab, nirsevimab | P03420 | 200 | 51 |
| SARS-CoV-2 Spike | Infectious | anti-SARS-CoV-2, anti-spike protein, bebtelovimab, sotrovimab, tixagevimab | P0DTC2 | 14 | 8 |
| Influenza HA | Infectious | anti-influenza hemagglutinin antibody, anti-HA broadly neutralizing | P03437 | 25 | 25 |
| HIV gp120 | Infectious | anti-HIV gp120, anti-HIV envelope, ibalizumab | P04578 | 7 | 2 |
| Ebola GP | Infectious | anti-Ebola glycoprotein, mAb114, atoltivimab, inmazeb | Q05320 | 7 | 2 |
| | | **Total** | | **7,672** | **3,829** |

**Script:** `scripts/fetch_patent_sequences.py`

---

## Source 3: ProteinBase Seed Data

**Database:** ProteinBase (literature-curated protein interaction pairs)

**Search strategy:** Pre-collected seed dataset of 300 antibody-antigen pairs
from published literature, curated from ProteinBase exports.

**Fields extracted:** binder name, binder sequence, binder type, antigen name,
antigen sequence, interaction label, affinity (KD), source reference (DOI/PMID)

**Script:** `scripts/build_antibody_antigen_candidates.py`

---

## Merge and Deduplication

All sources are merged by `scripts/merge_candidates.py`:

1. Records from each source are loaded and standardized to a common schema
2. Sequence pair keys (MD5 hash of binder + target sequences) are computed
3. Duplicate sequence pairs are removed (first occurrence kept)
4. Antigen names are normalized via `scripts/normalize_antigens.py`
5. Quality flags are computed (sequence contains X, invalid chars, etc.)
6. Binder types are classified and reviewed
7. Affinity values are normalized to nM (normalized_KD_nM)

**Output:** `antibody_antigen_master_v2.csv` (app-facing master table)

---

## PubMed Literature References

**Database:** PubMed / NCBI E-utilities

**Search strategy:** For each unique antigen in the master table, PubMed is
searched for relevant antibody engineering literature.

**Query template:** `{antigen_name} antibody`

**Parameters:** retmax=5 per target, max-targets=30

**Script:** `scripts/fetch_pubmed_references.py`

**Output:** `pubmed_references.csv` (120 rows covering 24 antigens)
