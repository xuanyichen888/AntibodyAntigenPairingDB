#!/usr/bin/env python3
"""
Step 1b: Fetch patent-deposited antibody sequences from NCBI GenBank PAT division.

Searches NCBI Entrez protein database for antibody variable regions targeting
specific antigens, then fetches canonical antigen sequences from UniProt.

Output: data/patent/patent_ab_ag_pairs.csv
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "data" / "patent"
OUTPUT_PATH = OUTPUT_DIR / "patent_ab_ag_pairs.csv"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"

TARGETS = {
    # --- Immune checkpoint ---
    "PD-1": {"query": "anti-PD-1 OR anti-PDCD1 OR nivolumab OR pembrolizumab OR cemiplimab", "uniprot": "Q15116"},
    "PD-L1": {"query": "anti-PD-L1 OR anti-CD274 OR atezolizumab OR durvalumab OR avelumab", "uniprot": "Q9NZQ7"},
    "CTLA-4": {"query": "anti-CTLA-4 OR anti-CTLA4 OR ipilimumab OR tremelimumab", "uniprot": "P16410"},
    "TIGIT": {"query": "anti-TIGIT OR tiragolumab OR vibostolimab OR domvanalimab", "uniprot": "Q495A1"},
    "LAG-3": {"query": "anti-LAG-3 OR anti-LAG3 OR relatlimab OR favezelimab", "uniprot": "P18627"},
    "TIM-3": {"query": "anti-TIM-3 OR anti-TIM3 OR anti-HAVCR2 OR cobolimab OR sabatolimab", "uniprot": "Q8TDQ0"},
    "VISTA": {"query": "anti-VISTA OR anti-VSIR OR anti-B7-H5", "uniprot": "Q9H7M9"},
    "BTLA": {"query": "anti-BTLA OR anti-CD272", "uniprot": "Q7Z6A9"},
    "CD47": {"query": "anti-CD47 OR magrolimab OR lemzoparlimab", "uniprot": "Q08722"},
    "SIRPa": {"query": "anti-SIRPa OR anti-SIRPA OR anti-CD172a", "uniprot": "P78324"},
    # --- Oncology surface targets ---
    "HER2": {"query": "anti-HER2 OR anti-ERBB2 OR trastuzumab OR pertuzumab OR margetuximab", "uniprot": "P04626"},
    "EGFR": {"query": "anti-EGFR OR anti-ERBB1 OR cetuximab OR panitumumab OR necitumumab", "uniprot": "P00533"},
    "HER3": {"query": "anti-HER3 OR anti-ERBB3 OR patritumab", "uniprot": "P21860"},
    "CD20": {"query": "anti-CD20 OR anti-MS4A1 OR rituximab OR obinutuzumab OR ofatumumab", "uniprot": "P11836"},
    "CD19": {"query": "anti-CD19 OR blinatumomab OR loncastuximab OR tafasitamab", "uniprot": "P15391"},
    "CD3": {"query": "anti-CD3 OR anti-CD3E OR blinatumomab OR teplizumab", "uniprot": "P07766"},
    "CD38": {"query": "anti-CD38 OR daratumumab OR isatuximab OR mezagitamab", "uniprot": "P28907"},
    "BCMA": {"query": "anti-BCMA OR anti-TNFRSF17 OR teclistamab OR elranatamab", "uniprot": "Q02223"},
    "GPRC5D": {"query": "anti-GPRC5D OR talquetamab", "uniprot": "Q9NZD1"},
    "FcRH5": {"query": "anti-FcRH5 OR anti-FCRL5 OR cevostamab", "uniprot": "Q96RD9"},
    "CD33": {"query": "anti-CD33 OR gemtuzumab", "uniprot": "P20138"},
    "CD22": {"query": "anti-CD22 OR inotuzumab OR moxetumomab", "uniprot": "P20273"},
    "CD30": {"query": "anti-CD30 OR anti-TNFRSF8 OR brentuximab", "uniprot": "P28908"},
    "TROP2": {"query": "anti-TROP2 OR anti-TACSTD2 OR sacituzumab OR datopotamab", "uniprot": "P09758"},
    "Nectin-4": {"query": "anti-Nectin-4 OR anti-PVRL4 OR enfortumab", "uniprot": "Q96NY8"},
    "CLDN18.2": {"query": "anti-claudin 18 OR anti-CLDN18 OR zolbetuximab", "uniprot": "P56856"},
    "DLL3": {"query": "anti-DLL3 OR rovalpituzumab OR tarlatamab", "uniprot": "Q9NYJ7"},
    "Mesothelin": {"query": "anti-mesothelin OR anti-MSLN OR anetumab", "uniprot": "Q13421"},
    "MET": {"query": "anti-MET OR anti-HGFR OR onartuzumab OR amivantamab", "uniprot": "P08581"},
    "PSMA": {"query": "anti-PSMA OR anti-FOLH1 OR capromab", "uniprot": "Q04609"},
    "FAP": {"query": "anti-FAP OR anti-fibroblast activation protein", "uniprot": "Q12884"},
    "GPC3": {"query": "anti-GPC3 OR anti-glypican-3 OR codrituzumab", "uniprot": "P51654"},
    "B7-H3": {"query": "anti-B7-H3 OR anti-CD276 OR enoblituzumab OR omburtamab", "uniprot": "Q5ZPR3"},
    "EpCAM": {"query": "anti-EpCAM OR anti-CD326 OR catumaxomab OR solitomab", "uniprot": "P16422"},
    # --- Cytokines & soluble targets ---
    "TNF-alpha": {"query": "anti-TNF OR anti-TNFA OR adalimumab OR infliximab OR golimumab OR certolizumab", "uniprot": "P01375"},
    "IL-6": {"query": "anti-IL-6 OR anti-IL6 OR siltuximab OR olokizumab", "uniprot": "P05231"},
    "IL-6R": {"query": "anti-IL-6R OR anti-IL6R OR tocilizumab OR sarilumab", "uniprot": "P08887"},
    "IL-17A": {"query": "anti-IL-17A OR anti-IL17A OR secukinumab OR ixekizumab OR bimekizumab", "uniprot": "Q16552"},
    "IL-17RA": {"query": "anti-IL-17RA OR anti-IL17RA OR brodalumab", "uniprot": "Q96F46"},
    "IL-4Ra": {"query": "anti-IL-4R OR anti-IL4R OR dupilumab", "uniprot": "P24394"},
    "IL-13": {"query": "anti-IL-13 OR anti-IL13 OR tralokinumab OR cendakimab", "uniprot": "P35225"},
    "IL-5": {"query": "anti-IL-5 OR anti-IL5 OR mepolizumab OR reslizumab", "uniprot": "P05113"},
    "IL-5Ra": {"query": "anti-IL-5R OR anti-IL5RA OR benralizumab", "uniprot": "Q01344"},
    "IL-23": {"query": "anti-IL-23 OR anti-IL23 OR guselkumab OR risankizumab OR tildrakizumab", "uniprot": "Q9NPF7"},
    "IL-12/23 p40": {"query": "anti-IL-12 OR anti-p40 OR ustekinumab", "uniprot": "P29460"},
    "IL-33": {"query": "anti-IL-33 OR anti-IL33 OR itepekimab OR tozorakimab", "uniprot": "O95760"},
    "IL-31": {"query": "anti-IL-31 OR anti-IL31 OR nemolizumab", "uniprot": "Q6EBC2"},
    "TSLP": {"query": "anti-TSLP OR tezepelumab", "uniprot": "Q969D9"},
    "IgE": {"query": "anti-IgE OR omalizumab OR ligelizumab", "uniprot": "P01854"},
    "VEGF": {"query": "anti-VEGF OR anti-VEGFA OR bevacizumab OR ranibizumab OR brolucizumab", "uniprot": "P15692"},
    "VEGFR2": {"query": "anti-VEGFR2 OR anti-KDR OR ramucirumab", "uniprot": "P35968"},
    "NGF": {"query": "anti-NGF OR tanezumab OR fasinumab OR fulranumab", "uniprot": "P01138"},
    "CGRP": {"query": "anti-CGRP OR fremanezumab OR galcanezumab OR eptinezumab", "uniprot": "P06881"},
    "CGRPr": {"query": "anti-CGRP receptor OR erenumab OR anti-CALCRL", "uniprot": "Q16602"},
    "PCSK9": {"query": "anti-PCSK9 OR alirocumab OR evolocumab", "uniprot": "Q8NBP7"},
    "C5": {"query": "anti-C5 OR eculizumab OR ravulizumab", "uniprot": "P01031"},
    "RANKL": {"query": "anti-RANKL OR anti-TNFSF11 OR denosumab", "uniprot": "O14788"},
    "Ang2": {"query": "anti-angiopoietin-2 OR anti-ANGPT2 OR faricimab", "uniprot": "O15123"},
    "BAFF": {"query": "anti-BAFF OR anti-TNFSF13B OR belimumab", "uniprot": "Q9Y275"},
    # --- Infectious disease ---
    "RSV F": {"query": "anti-RSV OR anti-respiratory syncytial virus OR palivizumab OR nirsevimab", "uniprot": "P03420"},
    "SARS-CoV-2 Spike": {"query": "anti-SARS-CoV-2 OR anti-spike protein OR bebtelovimab OR sotrovimab OR tixagevimab", "uniprot": "P0DTC2"},
    "Influenza HA": {"query": "anti-influenza hemagglutinin antibody OR anti-HA broadly neutralizing", "uniprot": "P03437"},
    "HIV gp120": {"query": "anti-HIV gp120 OR anti-HIV envelope OR ibalizumab", "uniprot": "P04578"},
    "Ebola GP": {"query": "anti-Ebola glycoprotein OR mAb114 OR atoltivimab OR inmazeb", "uniprot": "Q05320"},
}

FIELDNAMES = [
    "ncbi_accession", "binder_name", "binder_type", "binder_sequence",
    "antigen_name", "antigen_uniprot_id", "antigen_sequence",
    "patent_id", "source",
]

AA_CHARS = set("ACDEFGHIKLMNPQRSTVWYXBZU")


def parse_accession(header: str) -> str:
    token = header.split()[0]
    parts = token.split("|")
    for part in parts:
        if re.fullmatch(r"[A-Z]{2,4}_?\d+(?:\.\d+)?", part):
            return part
    for part in parts:
        if part and part.lower() not in {"gb", "emb", "dbj", "pat", "ref", "sp", "tr", "us"}:
            return part
    return token


def extract_patent_id(title: str) -> str:
    match = re.search(r"\b(?:US|WO|EP)\s*\d+[A-Z0-9 ]{0,8}", title, re.IGNORECASE)
    if match:
        return re.sub(r"\s+", " ", match.group(0)).strip(" .;,")
    match = re.search(r"Patent:\s*([^,\]\s]+)", title, re.IGNORECASE)
    return match.group(1) if match else ""


def classify_binder(title: str, seq: str = "") -> str:
    t = title.lower()
    if any(kw in t for kw in ("nanobody", "vhh", "single domain", "single-domain", "sdab")):
        return "nanobody"
    if any(kw in t for kw in ("scfv", "single chain", "single-chain")):
        return "scFv"
    if any(kw in t for kw in ("heavy chain", " vh ")):
        return "heavy chain"
    if any(kw in t for kw in ("light chain", " vl ")):
        return "light chain"
    if "fab" in t:
        return "Fab"
    if seq:
        slen = len(seq)
        if 200 <= slen <= 300:
            return "scFv"
        if 80 <= slen <= 150:
            return "antibody fragment"
    return "antibody"


def ncbi_search(query: str, retmax: int = 200) -> list[str]:
    full_query = f'({query}) AND {AB_FILTER} AND "PAT"[Division]'
    params = urllib.parse.urlencode({
        "db": "protein",
        "term": full_query,
        "retmax": retmax,
        "retmode": "json",
    })
    url = f"{ESEARCH_URL}?{params}"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                data = json.loads(resp.read())
            return data.get("esearchresult", {}).get("idlist", [])
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)
    return []


def ncbi_fetch_fasta(gis: list[str]) -> list[dict]:
    if not gis:
        return []
    entries = []
    batch_size = 100
    for i in range(0, len(gis), batch_size):
        batch = gis[i:i + batch_size]
        params = urllib.parse.urlencode({
            "db": "protein",
            "id": ",".join(batch),
            "rettype": "fasta",
            "retmode": "text",
        })
        url = f"{EFETCH_URL}?{params}"
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as resp:
                    text = resp.read().decode()
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2)

        for block in text.strip().split(">"):
            block = block.strip()
            if not block:
                continue
            lines = block.split("\n")
            header = lines[0]
            seq = "".join(line.strip() for line in lines[1:])
            accession = parse_accession(header)
            entries.append({"accession": accession, "title": header, "sequence": seq})
        if i + batch_size < len(gis):
            time.sleep(0.5)
    return entries


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    url = f"{UNIPROT_URL}/{uniprot_id}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            text = resp.read().decode()
        lines = text.strip().split("\n")
        return "".join(line.strip() for line in lines[1:])
    except Exception:
        return ""


AB_FILTER = (
    "(antibody OR immunoglobulin OR VH OR VHH OR VL OR scFv OR nanobody "
    "OR Fab OR heavy chain OR light chain OR CDR OR variable region)"
)


def is_antibody_sequence(seq: str, title: str) -> bool:
    if len(seq) < 80 or len(seq) > 600:
        return False
    return all(c in AA_CHARS for c in seq.upper())


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    antigen_seqs: dict[str, str] = {}
    print("Fetching antigen sequences from UniProt...")
    for name, info in TARGETS.items():
        seq = fetch_uniprot_sequence(info["uniprot"])
        if seq:
            antigen_seqs[name] = seq
            print(f"  {name}: {len(seq)} aa")
        else:
            print(f"  {name}: FAILED")
        time.sleep(0.3)

    rows: list[dict] = []
    for target_name, info in TARGETS.items():
        ag_seq = antigen_seqs.get(target_name, "")
        if not ag_seq:
            continue
        try:
            gis = ncbi_search(info["query"])
            print(f"  {target_name}: {len(gis)} hits")
        except Exception as e:
            print(f"  {target_name}: search failed - {e}")
            continue
        time.sleep(1.0)

        if not gis:
            continue
        try:
            entries = ncbi_fetch_fasta(gis)
        except Exception as e:
            print(f"  {target_name}: fetch failed - {e}")
            continue
        time.sleep(1.0)

        for entry in entries:
            if is_antibody_sequence(entry["sequence"], entry["title"]):
                patent_id = extract_patent_id(entry["title"])
                rows.append({
                    "ncbi_accession": entry["accession"],
                    "binder_name": entry["title"].split("|")[-1].strip()[:200] if "|" in entry["title"] else entry["title"][:200],
                    "binder_type": classify_binder(entry["title"], entry["sequence"]),
                    "binder_sequence": entry["sequence"],
                    "antigen_name": target_name,
                    "antigen_uniprot_id": info["uniprot"],
                    "antigen_sequence": ag_seq,
                    "patent_id": patent_id,
                    "source": f"GenBank:{entry['accession']}",
                })

    seen = set()
    deduped = []
    for r in rows:
        key = r["ncbi_accession"]
        if key not in seen:
            seen.add(key)
            deduped.append(r)
    rows = deduped

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {OUTPUT_PATH} ({len(rows)} rows, {len(seen)} unique accessions)")
    types = Counter(r["binder_type"] for r in rows)
    targets = Counter(r["antigen_name"] for r in rows)
    print(f"Binder types: {dict(types)}")
    print(f"Targets: {dict(targets)}")


if __name__ == "__main__":
    main()
