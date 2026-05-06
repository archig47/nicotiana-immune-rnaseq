import os

# File paths
hcrv_file = "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_HCRV_vs_CK.csv"
tswv_file = "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_TSWV_vs_CK.csv"
th_file   = "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_TH_vs_CK.csv"
promoter_fa = "/users/fyp/fyp5/project/PRJNA945175/promoters/all_PRR_promoters_2kb.fa"
outdir = "/users/fyp/fyp5/project/PRJNA945175/promoters/subsets"
os.makedirs(outdir, exist_ok=True)

def read_degs(filepath):
    up, down = set(), set()
    with open(filepath) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(',')
            gene = parts[6].strip().strip('"')
            lfc = float(parts[1])
            if lfc > 0:
                up.add(gene)
            else:
                down.add(gene)
    return up, down

hcrv_up, hcrv_down = read_degs(hcrv_file)
tswv_up, tswv_down = read_degs(tswv_file)
th_up,   th_down   = read_degs(th_file)

# Key subsets
subsets = {
    # Condition-specific upregulated
    "HCRV_up": hcrv_up,
    "TSWV_up": tswv_up,
    "TH_up":   th_up,
    "HCRV_down": hcrv_down,
    "TSWV_down": tswv_down,
    "TH_down":   th_down,

    # The 53 PRRs induced in single infections but suppressed in TH co-infection
    # = upregulated in HCRV or TSWV, but downregulated in TH
    "induced_single_suppressed_TH": (hcrv_up | tswv_up) & th_down,

    # TSWV-specific suppressors (down in TSWV, not in HCRV)
    "TSWV_specific_down": tswv_down - hcrv_down,

    # Universal responders across all 3 conditions
    "all3_up":   hcrv_up & tswv_up & th_up,
    "all3_down": hcrv_down & tswv_down & th_down,
}

# Print summary
for name, genes in subsets.items():
    print(f"{name}: {len(genes)} genes")

# Load promoter sequences
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

print(f"\nPromoters loaded: {len(promoters)}")

# Write subset FASTA files
for name, genes in subsets.items():
    if len(genes) < 5:
        print(f"  SKIPPING {name} — too few genes ({len(genes)})")
        continue
    outpath = os.path.join(outdir, f"{name}_promoters_2kb.fa")
    found = 0
    with open(outpath, 'w') as out:
        for gene in sorted(genes):
            if gene in promoters:
                out.write(f">{gene}\n{promoters[gene]}\n")
                found += 1
    print(f"  Written {found} -> {outpath}")
