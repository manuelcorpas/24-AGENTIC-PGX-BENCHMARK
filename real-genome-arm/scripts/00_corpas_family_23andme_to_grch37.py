#!/usr/bin/env python3
"""
European family arm — convert the public Corpas family 23andMe SNP-chip genotypes
(figshare: 23andMe SNP chip genotype data) into a single multi-sample GRCh37 VCF for
PyPGx. The 23andMe raw data is GRCh36/hg18, so each SNP is lifted hg18->hg19 and the
reference allele is set from the hg19 genome (so REF/ALT/genotype orientation is
correct and PyPGx matches).

Requires: pyliftover, pysam; the hg18ToHg19 chain (UCSC) and an hg19 reference fasta
(indexed). Member files are the per-person 23andMe exports (rsid, chromosome, position,
genotype), autosomes used.

Usage:
  00_corpas_family_23andme_to_grch37.py <chain.gz> <hg19.fa> <out.vcf> <member1.txt> [member2.txt ...]
"""
import sys
import pysam
from pyliftover import LiftOver

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AUT = set(str(i) for i in range(1, 23))

def main():
    chain, fasta, out = sys.argv[1], sys.argv[2], sys.argv[3]
    member_files = sys.argv[4:]
    members = [f.rsplit("/", 1)[-1].rsplit(".", 1)[0] for f in member_files]
    lo = LiftOver(chain)
    fa = pysam.FastaFile(fasta)
    fa_chr_prefix = fa.references[0].startswith("chr")

    data = {}  # (chrom, pos_hg18, rsid) -> {member: genotype}
    for f, m in zip(member_files, members):
        for line in open(f):
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[1] not in AUT:
                continue
            data.setdefault((p[1], int(p[2]), p[0]), {})[m] = p[3]

    recs = []
    for (chrom, pos36, rsid), gts in data.items():
        r = lo.convert_coordinate("chr" + chrom, pos36 - 1)
        if not r:
            continue
        c37 = r[0][0].replace("chr", "")
        pos37 = r[0][1] + 1
        if c37 != chrom:
            continue
        ref_contig = ("chr" + chrom) if fa_chr_prefix else chrom
        try:
            ref = fa.fetch(ref_contig, pos37 - 1, pos37).upper()
        except Exception:
            continue
        if ref not in "ACGT":
            continue
        alleles = set(); vg = {}
        for m, g in gts.items():
            if len(g) == 2 and all(ch in "ACGT" for ch in g):
                vg[m] = g; alleles |= set(g)
        if not vg:
            continue
        if ref not in alleles:  # try opposite strand
            ca = {COMP[a] for a in alleles}
            if ref in ca:
                vg = {m: "".join(COMP[ch] for ch in g) for m, g in vg.items()}
                alleles = ca
            else:
                continue
        nonref = alleles - {ref}
        if len(nonref) != 1:
            continue
        alt = nonref.pop()
        out_gt = {m: ("0/0" if g.count(alt) == 0 else ("0/1" if g.count(alt) == 1 else "1/1"))
                  for m, g in vg.items()}
        recs.append((int(chrom), pos37, rsid, ref, alt, out_gt))

    recs.sort()
    with open(out, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in range(1, 23):
            f.write(f"##contig=<ID={c}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(members) + "\n")
        for chrom, pos, rsid, ref, alt, out_gt in recs:
            cells = "\t".join(out_gt.get(m, "./.") for m in members)
            f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{cells}\n")
    print(f"wrote {len(recs)} variants for {len(members)} members -> {out}")
    print("then: bgzip + tabix index, and run 02_call_diplotypes_pypgx.sh")

if __name__ == "__main__":
    main()
