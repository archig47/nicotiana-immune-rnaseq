import re
import os

# Paths
genome    = "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.softmasked.genome.fasta"
gff       = "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3"
tracker   = "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.PRRtracker_classified.tsv"
outdir    = "/users/fyp/fyp5/project/PRJNA945175/promoters"
prr_all = set()
with open(tracker) as fh:
    header = next(fh)
    for line in fh:
        parts = line.strip().split()
        if parts:
            gene_id = parts[0].strip().replace('-mRNA', '')
            prr_all.add(gene_id)

print(f"PRR IDs loaded from tracker: {len(prr_all)}")
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
        m = re.search(r'ID=([^;]+)', attrs)
        if not m:
            continue
        gene_id = m.group(1)
        gene_coords[gene_id] = (chrom, start, end, strand)

print(f"Genes parsed from GFF: {len(gene_coords)}")

# Load genome
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

def extract_promoter(gene_id, upstream=2000):
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

# Write all PRR promoters
out_path = os.path.join(outdir, "all_PRR_promoters_2kb.fa")
found = 0
missing = 0
with open(out_path, 'w') as out:
    for gene_id in sorted(prr_all):
        seq = extract_promoter(gene_id)
        if seq:
            out.write(f">{gene_id}\n{seq}\n")
            found += 1
        else:
            missing += 1

print(f"PRR promoters written: {found} -> {out_path}")
print(f"Missing from GFF: {missing}")
