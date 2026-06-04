import re
import os

# Paths
genome = "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.softmasked.genome.fasta"
gff    = "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3"
outdir = "/users/fyp/fyp5/project/PRJNA945175/promoters"

# DEG NLR IDs across all three conditions
nlr_ids = set([
    "NbT2T02g01323","NbT2T03g03424","NbT2T04g02568","NbT2T04g02613",
    "NbT2T05g01385","NbT2T06g03156","NbT2T06g03387","NbT2T08g02564",
    "NbT2T09g00151","NbT2T10g00750","NbT2T13g01835","NbT2T14g00437",
    "NbT2T15g02674","NbT2T16g00401","NbT2T16g03031","NbT2T17g00315",
    "NbT2T17g02644","NbT2T18g00191","NbT2T19g01553","NbT2T03g00219",
    "NbT2T04g02572","NbT2T13g02081","NbT2T15g02837","NbT2T18g01444",
    "NbT2T19g01736","NbT2T02g02855","NbT2T05g01392","NbT2T07g02403",
    "NbT2T09g01711","NbT2T12g02670","NbT2T12g02679","NbT2T13g02756",
    "NbT2T15g01346","NbT2T17g00188"
])

# PRR DEG IDs — collect from all three condition files
prr_files = [
    "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_HCRV_vs_CK.csv",
    "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_TSWV_vs_CK.csv",
    "/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs32_PRR_TH_vs_CK.csv",
]

prr_ids = set()
for f in prr_files:
    with open(f) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.strip().split(',')
            gene = parts[6].strip().strip('"')
            prr_ids.add(gene)

print(f"NLR DEGs: {len(nlr_ids)}")
print(f"PRR DEGs: {len(prr_ids)}")

# Parse GFF3 for gene coordinates
gene_coords = {}
with open(gff) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'gene':
            continue
        chrom  = parts[0]
        start  = int(parts[3])
        end    = int(parts[4])
        strand = parts[6]
        attrs  = parts[8]
        # extract ID
        m = re.search(r'ID=([^;]+)', attrs)
        if not m:
            continue
        gene_id = m.group(1)
        gene_coords[gene_id] = (chrom, start, end, strand)

print(f"Genes parsed from GFF: {len(gene_coords)}")

# Load genome into dict
genome_seqs = {}
current_chrom = None
current_seq = []
with open(genome) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            if current_chrom:
                genome_seqs[current_chrom] = ''.join(current_seq)
            current_chrom = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_chrom:
        genome_seqs[current_chrom] = ''.join(current_seq)

print(f"Chromosomes loaded: {len(genome_seqs)}")

def rev_comp(seq):
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]

def extract_promoter(gene_id, gene_set, upstream=2000):
    if gene_id not in gene_set:
        return None
    if gene_id not in gene_coords:
        print(f"  WARNING: {gene_id} not in GFF")
        return None
    chrom, start, end, strand = gene_coords[gene_id]
    chrom_seq = genome_seqs.get(chrom)
    if chrom_seq is None:
        print(f"  WARNING: {chrom} not in genome")
        return None
    if strand == '+':
        prom_start = max(0, start - 1 - upstream)
        prom_end   = start - 1
        seq = chrom_seq[prom_start:prom_end]
    else:
        prom_start = end
        prom_end   = min(len(chrom_seq), end + upstream)
        seq = rev_comp(chrom_seq[prom_start:prom_end])
    return seq

# Write NLR promoters
nlr_out = os.path.join(outdir, "DEG_NLR_promoters_2kb.fa")
found_nlr = 0
with open(nlr_out, 'w') as out:
    for gene_id in sorted(nlr_ids):
        seq = extract_promoter(gene_id, nlr_ids)
        if seq:
            out.write(f">{gene_id}\n{seq}\n")
            found_nlr += 1
print(f"NLR promoters written: {found_nlr} -> {nlr_out}")

# Write PRR promoters
prr_out = os.path.join(outdir, "DEG_PRR_promoters_2kb.fa")
found_prr = 0
with open(prr_out, 'w') as out:
    for gene_id in sorted(prr_ids):
        seq = extract_promoter(gene_id, prr_ids)
        if seq:
            out.write(f">{gene_id}\n{seq}\n")
            found_prr += 1
print(f"PRR promoters written: {found_prr} -> {prr_out}")
