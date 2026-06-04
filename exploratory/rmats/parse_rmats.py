#!/usr/bin/env python3
"""
Parse rMATS results for NLR and PRR genes.
Filters for FDR < 0.05 and |ΔPSI| > 0.1 (JC counts).
Outputs per-gene summary and per-event detail tables.
"""
import os, sys
import pandas as pd
import numpy as np

BASE    = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
RMATS   = "/users/fyp/fyp5/project/PRJNA945175/rmats_v2"
NLR_IDS = os.path.join(BASE, "nlr_prr/nlr_gene_ids.txt")
PRR_XLS = os.path.join(BASE, "nlr_prr/nlr_prr_full_landscape.xlsx")
NLR_CLS = os.path.join(BASE, "nlr_prr/NLR_phylogeneti_classification.xlsx")
OUT_DIR = os.path.join(BASE, "figures")

os.makedirs(OUT_DIR, exist_ok=True)

# ── Gene lists ─────────────────────────────────────────────────────────────
nlr_ids = set(open(NLR_IDS).read().strip().split())

prr_df  = pd.read_excel(PRR_XLS, sheet_name="PRR_landscape")
prr_ids = set(prr_df["gene_id"].dropna().astype(str))

print(f"NLR genes: {len(nlr_ids)}")
print(f"PRR genes: {len(prr_ids)}")

# NLR subclass (for annotation)
try:
    nlr_cls_df = pd.read_excel(NLR_CLS)
    # Find gene_id and subclass columns
    gene_col = [c for c in nlr_cls_df.columns if "gene" in c.lower() or "id" in c.lower()][0]
    cls_col  = [c for c in nlr_cls_df.columns if "class" in c.lower() or "subclass" in c.lower() or "clade" in c.lower()][0]
    # Strip -mRNA suffix from IDs in classification file
    nlr_subclass = dict(zip(
        nlr_cls_df[gene_col].astype(str).str.replace(r'-mRNA.*$', '', regex=True),
        nlr_cls_df[cls_col].astype(str)
    ))
    print(f"NLR subclass map: {len(nlr_subclass)} entries, col used: {gene_col} / {cls_col}")
except Exception as e:
    print(f"WARNING: NLR subclass load failed: {e}")
    nlr_subclass = {}

# PRR ecto classification (same logic as figures)
ACCESSORY = {"EGF","PAN","SDOM","NEW","RCC1","TNFR","GP_PDE","UNIDENTIFIED"}

def get_ecto(rules):
    if not isinstance(rules, str) or pd.isna(rules):
        return "Other"
    tokens = rules.split(",")[1:]
    primary = next((t for t in tokens if t not in ACCESSORY), None)
    if primary is None:         return "Other"
    if primary == "LRR":        return "LRR"
    if primary == "LysM":       return "LysM"
    if primary == "WAK":        return "WAK"
    if primary in ("MAL","SPARK"):            return "Malectin/CrRLK"
    if primary in ("BLEC","LLEC","CLEC","GLEC","GNK2"): return "Lectin"
    return "Other"

prr_ecto = {row["gene_id"]: get_ecto(row.get("matched_rules","")) for _, row in prr_df.iterrows()}

COMPARISONS = [
    "HCRV_1dpi","HCRV_7dpi","HCRV_14dpi",
    "TSWV_1dpi","TSWV_7dpi","TSWV_14dpi",
    "TH_1dpi",  "TH_7dpi",  "TH_14dpi"
]

EVENT_TYPES = ["SE","A5SS","A3SS","MXE","RI"]

FDR_THRESH  = 0.05
DPSI_THRESH = 0.1

# ── Parse ──────────────────────────────────────────────────────────────────
records = []
for comp in COMPARISONS:
    for evt in EVENT_TYPES:
        fpath = os.path.join(RMATS, comp, f"{evt}.MATS.JC.txt")
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath, sep="\t")
        if df.empty:
            continue
        # GeneID column is the gene name field in rMATS output
        df = df.rename(columns={"GeneID": "gene_id"})
        # rMATS gene IDs from GTF are "Nb01g..." — convert to NbT2T01g... format
        df["gene_id"] = (df["gene_id"].astype(str)
                         .str.strip('"')
                         .str.replace(r'^Nb(\d)', r'NbT2T\1', regex=True))
        # Filter for NLR or PRR
        is_nlr = df["gene_id"].isin(nlr_ids)
        is_prr = df["gene_id"].isin(prr_ids)
        sub = df[is_nlr | is_prr].copy()
        if sub.empty:
            continue
        sub["gene_class"] = "PRR"
        sub.loc[is_nlr[sub.index], "gene_class"] = "NLR"
        sub["comparison"] = comp
        sub["event_type"] = evt
        sub["FDR"]  = pd.to_numeric(sub["FDR"],  errors="coerce")
        sub["DPSI"] = pd.to_numeric(sub["IncLevelDifference"], errors="coerce")
        # Significant: FDR<0.05 AND |ΔPSI|>0.1
        sub["significant"] = (sub["FDR"] < FDR_THRESH) & (sub["DPSI"].abs() > DPSI_THRESH)
        records.append(sub)

if not records:
    print("ERROR: no records found — check gene ID format in SE files")
    sys.exit(1)

all_df = pd.concat(records, ignore_index=True)
print(f"\nTotal rMATS rows in NLR/PRR genes: {len(all_df)}")
print(f"Significant (FDR<{FDR_THRESH}, |ΔPSI|>{DPSI_THRESH}): {all_df['significant'].sum()}")

sig_df = all_df[all_df["significant"]].copy()

# Add annotation
sig_df["subclass_ecto"] = sig_df.apply(
    lambda r: nlr_subclass.get(r["gene_id"], "Unknown") if r["gene_class"] == "NLR"
              else prr_ecto.get(r["gene_id"], "Other"),
    axis=1
)

# ── Summary tables ─────────────────────────────────────────────────────────
# 1. Per comparison × gene_class × event_type count
summary = (sig_df.groupby(["comparison","gene_class","event_type"])
           .size().reset_index(name="n_sig_events"))
print("\n=== Significant AS events per comparison/class/event_type ===")
print(summary.to_string(index=False))

# 2. Unique significant genes per comparison × gene_class
unique_genes = (sig_df.groupby(["comparison","gene_class"])["gene_id"]
                .nunique().reset_index(name="n_genes_with_sig_AS"))
print("\n=== Unique genes with ≥1 sig AS event ===")
print(unique_genes.to_string(index=False))

# 3. Genes significant in multiple comparisons (most frequent)
gene_comp_count = (sig_df.groupby(["gene_class","gene_id"])["comparison"]
                   .nunique().reset_index(name="n_comparisons"))
gene_comp_count["subclass_ecto"] = gene_comp_count.apply(
    lambda r: nlr_subclass.get(r["gene_id"], "Unknown") if r["gene_class"] == "NLR"
              else prr_ecto.get(r["gene_id"], "Other"),
    axis=1
)
top = gene_comp_count.sort_values(["gene_class","n_comparisons"], ascending=[True,False])
print("\n=== Genes with sig AS in most comparisons (top 20) ===")
print(top.head(20).to_string(index=False))

# ── Save detail tables ────────────────────────────────────────────────────
keep_cols = ["comparison","gene_class","event_type","gene_id","subclass_ecto",
             "chr","strand","FDR","DPSI","IncLevel1","IncLevel2","significant"]
keep_cols = [c for c in keep_cols if c in sig_df.columns]

out_path = os.path.join(OUT_DIR, "rmats_nlr_prr_sig_events.csv")
sig_df[keep_cols].to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")

out_summary = os.path.join(OUT_DIR, "rmats_summary.csv")
unique_genes.to_csv(out_summary, index=False)
print(f"Saved: {out_summary}")

# ── SE (skipped exon) focus — most abundant event type ────────────────────
se_sig = sig_df[sig_df["event_type"] == "SE"]
print(f"\n=== SE events only: {len(se_sig)} significant ===")
print(se_sig.groupby(["comparison","gene_class"]).size().unstack(fill_value=0).to_string())
