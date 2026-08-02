#!/usr/bin/env python3
# aggregate_results.py — collect every incR2_<thresh>.tsv (+ _diff) into two tidy
# master tables for the boss: per-PRS incremental R2, and each PRS minus FELIX.
#   python3 aggregate_results.py <results_dir> <thresh>   # thresh e.g. 5e-8
import sys, os, glob
res, thr = sys.argv[1], sys.argv[2]
out1 = os.path.join(res, f"summary_{thr}.tsv")
out2 = os.path.join(res, f"diff_vs_felix_{thr}.tsv")

with open(out1, "w") as o1:
    o1.write("trait\tancestry\tscore\tpredictor\tldref\tN\tincR2\tCI_low\tCI_high\n")
    for f in sorted(glob.glob(os.path.join(res, "*", "*", f"incR2_{thr}.tsv"))):
        p = f.split(os.sep); trait, anc = p[-3], p[-2]
        with open(f) as fh:
            idx = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
            for line in fh:
                v = line.rstrip("\n").split("\t"); sc = v[idx["score"]]
                pred = sc.split("__")[0]; ld = sc.split("__")[1] if "__" in sc else ""
                o1.write("\t".join([trait, anc, sc, pred, ld, v[idx["N"]],
                                    v[idx["estimate"]], v[idx["CI_low"]], v[idx["CI_high"]]]) + "\n")

with open(out2, "w") as o2:
    o2.write("trait\tancestry\tcomparison\td_estimate\tCI_low\tCI_high\tp_boot\n")
    for f in sorted(glob.glob(os.path.join(res, "*", "*", f"incR2_{thr}_diff.tsv"))):
        p = f.split(os.sep); trait, anc = p[-3], p[-2]
        with open(f) as fh:
            idx = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
            for line in fh:
                v = line.rstrip("\n").split("\t")
                o2.write("\t".join([trait, anc] + [v[idx[c]] for c in
                         ["comparison", "d_estimate", "CI_low", "CI_high", "p_boot"]]) + "\n")
print("wrote", out1, "and", out2)
