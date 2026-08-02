#!/usr/bin/env python3
# build_geno2eid.py — write the id map bootstrap_r2.R expects: col1 = eid (matches
# --pheno_id and --covar_id 's'), col2 = genotype/sscore IID. If the sscore IIDs are
# already eids, write an EMPTY file (identity) so the caller omits --map.
#   build_geno2eid.py <sample_sscore> <pheno_tsv> <sample_id_mapping> <out_map>
import sys
sscore, pheno, mapping, out = sys.argv[1:5]

# genotype/sscore IIDs
G = set()
with open(sscore) as f:
    f.readline()
    for line in f:
        w = line.split()
        if w:
            G.add(w[0])

# phenotype eids (column 'eid')
P = set()
with open(pheno) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    j = hdr.index("eid") if "eid" in hdr else 0
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) > j:
            P.add(p[j])

if G and len(G & P) / len(G) > 0.5:
    open(out, "w").close()                       # identity -> caller omits --map
    print(f"IDENTITY: sscore IID already == eid ({len(G & P)}/{len(G)} overlap); no --map needed")
    sys.exit(0)

c1, c2 = [], []
with open(mapping) as f:
    for line in f:
        p = line.split()
        if len(p) >= 2:
            c1.append(p[0]); c2.append(p[1])
s1, s2 = set(c1), set(c2)
# the genotype column is whichever overlaps the sscore IIDs; the other is the eid column
if len(G & s2) >= len(G & s1):
    eidc, geno, og, oe = c1, c2, len(G & s2), len(P & s1)
else:
    eidc, geno, og, oe = c2, c1, len(G & s1), len(P & s2)
with open(out, "w") as o:
    for e, g in zip(eidc, geno):
        o.write(f"{e}\t{g}\n")                    # col1=eid, col2=genoIID
print(f"map -> {out}: genotype-col overlap {og}/{len(G)}, eid-col-vs-pheno {oe}/{len(P)}")
if og == 0:
    sys.stderr.write("WARN: genotype IIDs not found in either sample_id_mapping column!\n")
