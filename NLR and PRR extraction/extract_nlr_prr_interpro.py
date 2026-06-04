import json
import pandas as pd
from collections import defaultdict

INTERPRO_JSON = "/users/fyp/fyp5/project/genome/interpro_full/interpro/NbT2T.final_v12_proteins.fa.json"
NLR_TRACKER   = "/users/fyp/fyp5/project/genome/interpro_full/classification/NbT2T.final_v12_proteins.fa_NLRtracker_classified.tsv"
PRR_TRACKER   = "/users/fyp/fyp5/project/genome/interpro_full/classification/NbT2T.final_v12_proteins.fa_PRRtracker_classified.tsv"
DEG_XLSX      = "/users/fyp/fyp5/project/PRJNA945175/deseq2/nlr_prr_degs_timepoint.xlsx"
OUT_XLSX      = "/users/fyp/fyp5/project/PRJNA945175/deseq2/nlr_prr_full_landscape.xlsx"

print("Loading InterPro JSON...")
with open(INTERPRO_JSON) as f:
    data = json.load(f)
results = data["results"]
print(f"  {len(results)} proteins in JSON")

interpro_domains = defaultdict(list)
for entry in results:
    xref = entry.get("xref", [])
    if not xref:
        continue
    gene_id = xref[0]["id"].replace("-mRNA", "")
    for match in entry.get("matches", []):
        sig = match.get("signature", {})
        ipr_entry = sig.get("entry")
        acc  = sig.get("accession", "")
        name = sig.get("name", "")
        etype = sig.get("type", "")
        if ipr_entry:
            acc   = ipr_entry.get("accession", acc)
            name  = ipr_entry.get("name", name)
            etype = ipr_entry.get("type", etype)
        if acc:
            interpro_domains[gene_id].append({"accession": acc, "name": name, "type": etype})
print(f"  {len(interpro_domains)} genes with domain annotations")

def collapse_domains(gene_id):
    domains = interpro_domains.get(gene_id, [])
    if not domains:
        return pd.Series({"interpro_accessions": "", "interpro_names": "", "interpro_types": ""})
    seen = {}
    for d in domains:
        if d["accession"] not in seen:
            seen[d["accession"]] = d
    unique = list(seen.values())
    return pd.Series({
        "interpro_accessions": "; ".join(d["accession"] for d in unique),
        "interpro_names":      "; ".join(d["name"] or "" for d in unique),
        "interpro_types":      "; ".join(d["type"]       for d in unique),
    })

print("Loading tracker files...")
nlr = pd.read_csv(NLR_TRACKER, sep="\t")
prr = pd.read_csv(PRR_TRACKER, sep="\t")
nlr["gene_id"] = nlr["protein_id"].str.replace("-mRNA", "", regex=False)
prr["gene_id"] = prr["protein_id"].str.replace("-mRNA", "", regex=False)
print(f"  NLRtracker: {len(nlr)} genes")
print(f"  PRRtracker: {len(prr)} genes")

print("Loading DEG data...")
nlr_degs = pd.read_excel(DEG_XLSX, sheet_name="NLR_DEGs")
prr_degs = pd.read_excel(DEG_XLSX, sheet_name="PRR_DEGs")
deg_cols = ["gene_id", "mean_baseMean", "mean_abs_LFC", "n_sig", "direction", "peak_contrast"]
nlr_deg_info = nlr_degs[deg_cols].copy()
nlr_deg_info["is_DEG"] = True
prr_deg_info = prr_degs[deg_cols].copy()
prr_deg_info["is_DEG"] = True

print("Building NLR landscape...")
nlr_domains = nlr["gene_id"].apply(collapse_domains)
nlr_out = pd.concat([nlr, nlr_domains], axis=1)
nlr_out = nlr_out.merge(nlr_deg_info, on="gene_id", how="left")
nlr_out["is_DEG"] = nlr_out["is_DEG"].fillna(False)
tracker_nlr_cols = [c for c in nlr.columns if c != "gene_id"]
nlr_out = nlr_out[["gene_id"] + tracker_nlr_cols + ["interpro_accessions", "interpro_names", "interpro_types"] + ["is_DEG", "direction", "n_sig", "mean_baseMean", "mean_abs_LFC", "peak_contrast"]]

print("Building PRR landscape...")
prr_domains = prr["gene_id"].apply(collapse_domains)
prr_out = pd.concat([prr, prr_domains], axis=1)
prr_out = prr_out.merge(prr_deg_info, on="gene_id", how="left")
prr_out["is_DEG"] = prr_out["is_DEG"].fillna(False)
tracker_prr_cols = [c for c in prr.columns if c != "gene_id"]
prr_out = prr_out[["gene_id"] + tracker_prr_cols + ["interpro_accessions", "interpro_names", "interpro_types"] + ["is_DEG", "direction", "n_sig", "mean_baseMean", "mean_abs_LFC", "peak_contrast"]]

print("\n── NLR summary ──")
print(nlr_out["type"].value_counts().to_string())
print(f"\nNLRs with InterPro domains: {(nlr_out['interpro_accessions'] != '').sum()} / {len(nlr_out)}")
print(f"NLR DEGs: {nlr_out['is_DEG'].sum()} / {len(nlr_out)}")
print("\n── PRR summary ──")
print(prr_out["type"].value_counts().to_string())
print(f"\nPRRs with InterPro domains: {(prr_out['interpro_accessions'] != '').sum()} / {len(prr_out)}")
print(f"PRR DEGs: {prr_out['is_DEG'].sum()} / {len(prr_out)}")

print(f"\nWriting {OUT_XLSX}...")
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    nlr_out.to_excel(writer, sheet_name="NLR_landscape", index=False)
    prr_out.to_excel(writer, sheet_name="PRR_landscape", index=False)
print("Done.")
