import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(
    page_title="Antibody-Antigen PairingDB",
    page_icon="🧬",
    layout="wide",
)

PROJECT_ROOT = Path(__file__).parent
DATA_PATH = PROJECT_ROOT / "antibody_antigen_master_v2.csv"
LITERATURE_PATH = PROJECT_ROOT / "pubmed_references.csv"
V1_PATH = PROJECT_ROOT / "outputs" / "candidates_v1.csv"
V2_PATH = PROJECT_ROOT / "outputs" / "candidates_v2.csv"


@st.cache_data
def load_data(mtime: float):
    df = pd.read_csv(DATA_PATH, dtype=str).fillna("")
    df["seq_len_binder"] = df["binder_sequence"].str.len()
    df["seq_len_target"] = df["target_sequence"].str.len()
    return df


@st.cache_data
def load_literature(mtime: float):
    return pd.read_csv(LITERATURE_PATH, dtype=str).fillna("")


def search_df(df: pd.DataFrame, query: str) -> pd.DataFrame:
    if not query.strip():
        return df
    terms = query.lower().split()
    searchable = (
        df["binder_name"].str.lower() + " " +
        df["antigen_name"].str.lower() + " " +
        df["target_name"].str.lower() + " " +
        df["binder_type"].str.lower() + " " +
        df["source_type"].str.lower() + " " +
        df["evidence_text"].str.lower() + " " +
        df["source_reference"].str.lower()
    )
    mask = pd.Series(True, index=df.index)
    for term in terms:
        mask &= searchable.str.contains(term, na=False, regex=False)
    return df[mask]


st.title("Antibody-Antigen PairingDB")

if not DATA_PATH.exists():
    st.error(f"Data file not found: {DATA_PATH}")
    st.info("Run `python scripts/merge_candidates.py` first to generate the data.")
    st.stop()

mtime = DATA_PATH.stat().st_mtime
df = load_data(mtime)

# ── Sidebar filters ──────────────────────────────────────────────────────────
st.sidebar.header("Filters")

query = st.sidebar.text_input("Search (binder, antigen, source, evidence)", "")

source_options = sorted(df["source_type"].unique())
selected_sources = st.sidebar.multiselect("Source", source_options, default=source_options)

type_options = sorted(df["binder_type"].unique())
selected_types = st.sidebar.multiselect("Binder Type", type_options, default=type_options)

label_options = sorted(df["interaction_label"].unique())
selected_labels = st.sidebar.multiselect("Interaction Label", label_options, default=label_options)

conf_options = sorted(df["confidence"].unique())
selected_conf = st.sidebar.multiselect("Confidence", conf_options, default=conf_options)

antigen_list = sorted(df["antigen_name"].unique())
selected_antigens = st.sidebar.multiselect(
    f"Antigen ({len(antigen_list)} unique)",
    antigen_list,
    default=[],
    help="Leave empty to show all antigens",
)

min_binder_len = st.sidebar.slider("Min binder sequence length", 0, 600, 0)
max_binder_len = st.sidebar.slider("Max binder sequence length", 0, 2000, 2000)

# ── Apply filters ─────────────────────────────────────────────────────────────
filtered = df.copy()
filtered = search_df(filtered, query)
filtered = filtered[filtered["source_type"].isin(selected_sources)]
filtered = filtered[filtered["binder_type"].isin(selected_types)]
filtered = filtered[filtered["interaction_label"].isin(selected_labels)]
filtered = filtered[filtered["confidence"].isin(selected_conf)]

if selected_antigens:
    filtered = filtered[filtered["antigen_name"].isin(selected_antigens)]

filtered = filtered[
    (filtered["seq_len_binder"] >= min_binder_len) &
    (filtered["seq_len_binder"] <= max_binder_len)
]

# ── Summary metrics ───────────────────────────────────────────────────────────
col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Total Pairs", f"{len(filtered):,}")
col2.metric("Unique Antigens", filtered["antigen_name"].nunique())
col3.metric("Sources", filtered["source_type"].nunique())
col4.metric("Binder Types", filtered["binder_type"].nunique())
col5.metric("With Affinity", (filtered["affinity_value"] != "").sum())

# ── Distribution charts ──────────────────────────────────────────────────────
tab_table, tab_charts, tab_detail, tab_literature = st.tabs(
    ["Data Table", "Charts", "Pair Detail", "Literature"]
)

with tab_table:
    display_cols = [
        "record_id", "binder_name", "binder_type", "antigen_name",
        "interaction_label", "confidence", "source_type",
        "affinity_value", "affinity_unit",
        "seq_len_binder", "seq_len_target",
        "source_reference", "curation_status",
    ]
    st.dataframe(
        filtered[display_cols],
        use_container_width=True,
        height=600,
    )

    csv_data = filtered.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download filtered data (CSV)",
        csv_data,
        "antibody_antigen_filtered.csv",
        "text/csv",
    )

with tab_charts:
    chart_col1, chart_col2 = st.columns(2)

    with chart_col1:
        st.subheader("By Source")
        source_counts = filtered["source_type"].value_counts()
        st.bar_chart(source_counts)

        st.subheader("By Binder Type")
        type_counts = filtered["binder_type"].value_counts()
        st.bar_chart(type_counts)

    with chart_col2:
        st.subheader("Top 15 Antigens")
        ag_counts = filtered["antigen_name"].value_counts().head(15)
        st.bar_chart(ag_counts)

        st.subheader("Binder Sequence Length Distribution")
        len_data = filtered["seq_len_binder"].astype(int)
        st.bar_chart(len_data.value_counts().sort_index())

with tab_detail:
    st.subheader("View Pair Details")
    if len(filtered) == 0:
        st.info("No rows match the current filters.")
    else:
        pair_options = filtered.apply(
            lambda r: f"{r['record_id']} | {r['binder_name'][:40]} → {r['antigen_name'][:30]}",
            axis=1,
        ).tolist()
        selected_pair = st.selectbox("Select a pair", pair_options)
        if selected_pair:
            rid = selected_pair.split(" | ")[0]
            row = filtered[filtered["record_id"] == rid].iloc[0]

            detail_col1, detail_col2 = st.columns(2)
            with detail_col1:
                st.markdown("**Binder**")
                st.write(f"**Name:** {row['binder_name']}")
                st.write(f"**Type:** {row['binder_type']}")
                st.write(f"**Length:** {row['seq_len_binder']} aa")
                st.text_area("Binder Sequence", row["binder_sequence"], height=120, disabled=True)

            with detail_col2:
                st.markdown("**Antigen / Target**")
                st.write(f"**Antigen:** {row['antigen_name']}")
                st.write(f"**Target:** {row['target_name']}")
                st.write(f"**Length:** {row['seq_len_target']} aa")
                st.text_area("Target Sequence", row["target_sequence"], height=120, disabled=True)

            st.markdown("---")
            meta_col1, meta_col2, meta_col3 = st.columns(3)
            with meta_col1:
                st.write(f"**Interaction:** {row['interaction_label']}")
                st.write(f"**Confidence:** {row['confidence']}")
                st.write(f"**Curation:** {row['curation_status']}")
            with meta_col2:
                st.write(f"**Source:** {row['source_type']}")
                st.write(f"**Dataset:** {row.get('source_dataset', '')}")
                st.write(f"**Reference:** {row['source_reference']}")
            with meta_col3:
                if row["affinity_value"]:
                    st.write(f"**Affinity:** {row['affinity_value']} {row['affinity_unit']}")
                if row.get("evidence_text"):
                    st.write(f"**Evidence:** {row['evidence_text'][:200]}")

with tab_literature:
    st.subheader("PubMed References")
    if not LITERATURE_PATH.exists():
        st.info("Run `python3 scripts/fetch_pubmed_references.py` to generate PubMed references.")
    else:
        lit_mtime = LITERATURE_PATH.stat().st_mtime
        lit_df = load_literature(lit_mtime)
        visible_antigens = set(filtered["antigen_name"].unique())
        lit_filtered = lit_df[lit_df["antigen_name"].isin(visible_antigens)].copy()

        lit_query = st.text_input("Search literature", "", key="literature_search")
        if lit_query.strip():
            terms = lit_query.lower().split()
            searchable = (
                lit_filtered["antigen_name"].str.lower() + " " +
                lit_filtered["title"].str.lower() + " " +
                lit_filtered["journal"].str.lower() + " " +
                lit_filtered["authors"].str.lower() + " " +
                lit_filtered["doi"].str.lower() + " " +
                lit_filtered["pmid"].str.lower()
            )
            mask = pd.Series(True, index=lit_filtered.index)
            for term in terms:
                mask &= searchable.str.contains(term, na=False, regex=False)
            lit_filtered = lit_filtered[mask]

        lit_col1, lit_col2, lit_col3 = st.columns(3)
        lit_col1.metric("References", f"{len(lit_filtered):,}")
        lit_col2.metric("Antigens", lit_filtered["antigen_name"].nunique())
        lit_col3.metric("Journals", lit_filtered["journal"].nunique())

        lit_cols = [
            "antigen_name", "pmid", "title", "journal", "pubdate",
            "doi", "authors", "url",
        ]
        st.dataframe(lit_filtered[lit_cols], use_container_width=True, height=520)
        lit_csv = lit_filtered.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download literature data (CSV)",
            lit_csv,
            "pubmed_references_filtered.csv",
            "text/csv",
        )

# ── Footer ────────────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    f"**{len(df):,}** total pairs | "
    f"**{df['antigen_name'].nunique()}** antigens | "
    f"**{df['source_type'].nunique()}** sources"
)
if LITERATURE_PATH.exists():
    lit_count = len(pd.read_csv(LITERATURE_PATH, dtype=str))
    st.sidebar.markdown(f"**{lit_count:,}** PubMed references")

versions = []
if V1_PATH.exists():
    v1 = pd.read_csv(V1_PATH, dtype=str)
    versions.append(f"v1: {len(v1)} rows (seed)")
if V2_PATH.exists():
    v2 = pd.read_csv(V2_PATH, dtype=str)
    versions.append(f"v2: {len(v2)} rows (expanded)")
if versions:
    st.sidebar.markdown("**Versions:** " + " | ".join(versions))
