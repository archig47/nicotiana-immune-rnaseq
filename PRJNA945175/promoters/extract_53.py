import os

hcrv_deg = "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_HCRV_vs_CK.csv"
tswv_deg = "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_TSWV_vs_CK.csv"
th_all   = "/users/fyp/fyp5/project/PRJNA945175/deseq2/all_results32_TH_vs_CK.csv"
promoter_fa = "/users/fyp/fyp5/project/PRJNA945175/promoters/all_PRR_promoters_2kb.fa"
outdir = "/users/fyp/fyp5/project/PRJNA945175/promoters/subsets"
os.makedirs(outdir, exist_ok=True)

# Upregulated PRR DEGs in HCRV or TSWV
single_up = set()
for f in [hcrv_deg, tswv_deg]:
    with open(f) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(',')
            gene = parts[6].strip().strip('"')
            lfc = float(parts[1])
            if lfc > 0:
                single_up.add(gene)

print(f"PRRs upregulated in HCRV or TSWV: {len(single_up)}")

# All TH LFCs regardless of significance
th_lfc = {}
with open(th_all) as fh:
    next(fh)
    for line in fh:
        parts = line.strip().split(',')
        gene = parts[6].strip().strip('"')
        try:
            lfc = float(parts[1])
            th_lfc[gene] = lfc
        except:
            pass

# Find genes upregulated in single infections but negative LFC in TH
suppressed_in_TH = set()
for gene in single_up:
    if gene in th_lfc and th_lfc[gene] < 0:
        suppressed_in_TH.add(gene)

print(f"Of these, suppressed (LFC < 0) in TH: {len(suppressed_in_TH)}")

# Load promoters
promoters = {}
current_id = None
current_seq = []
with open(promoter_fa) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                promoters[current_id] = ''.join(current_seq)
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        promoters[current_id] = ''.join(current_seq)

# Write output
outpath = os.path.join(outdir, "induced_suppressed_TH_promoters_2kb.fa")
found = 0
with open(outpath, 'w') as out:
    for gene in sorted(suppressed_in_TH):
        if gene in promoters:
            out.write(f">{gene}\n{promoters[gene]}\n")
            found += 1

print(f"Written {found} -> {outpath}")
