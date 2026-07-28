#!/usr/bin/env python3
# Supplementary Table 9: European PRS prediction R2, FELIXassoc-local vs
# All-by-All-global per trait, with analytic 95% CIs. Writes .tsv (+ .xlsx).
import os, csv

PRS_DIR = os.path.dirname(os.path.abspath(__file__))
PREFER_ADJUSTED = False
RAW = "results/r2_4way.tsv"
ADJ = "results/r2_4way_adjusted_EUR.tsv"

TRAITS = [  # (folder id, display name)
    ("3006923", "Alanine aminotransferase (ALT)"),
    ("3007070", "HDL cholesterol"),
    ("3009744", "MCHC"),
    ("3013721", "Aspartate aminotransferase (AST)"),
    ("3022192", "Triglycerides"),
    ("3024929", "Platelet count"),
    ("3027114", "Total cholesterol"),
    ("3028288", "LDL cholesterol"),
    ("3035995", "Alkaline phosphatase (ALP)"),
    ("BMI",     "BMI"),
    ("height",  "Height"),
]

def read_r2(idd):
    fa = os.path.join(PRS_DIR, idd, ADJ)
    fr = os.path.join(PRS_DIR, idd, RAW)
    f = fa if (PREFER_ADJUSTED and os.path.exists(fa)) else fr
    g = l = n = None
    with open(f) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0] == "":            # header line
                continue
            if row[0] == "EUR-global": g, n = float(row[2]), float(row[1])
            if row[0] == "EUR-local":  l, n = float(row[2]), float(row[1])
    return g, l, n

def read_training_N():
    """AoU EUR local training N from phenotype_local_sample_sizes.tsv (col anc3/EUR)."""
    p = os.path.join(PRS_DIR, "phenotype_local_sample_sizes.tsv")
    out = {}
    if not os.path.exists(p):
        return out
    with open(p) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr = rows[0]
    # find an EUR-local column heuristically
    eur_col = next((i for i, c in enumerate(hdr) if "EUR" in c.upper()), None)
    id_col = 0
    for r in rows[1:]:
        if len(r) > (eur_col or 0):
            out[r[id_col]] = r[eur_col] if eur_col is not None else ""
    return out

def se_r2(r2, n):
    """large-sample SE of R2 (squared correlation): 2*sqrt(R2)*(1-R2)/sqrt(N)."""
    return 2.0 * (max(r2, 0.0) ** 0.5) * (1.0 - r2) / (n ** 0.5)

def ci(r2, n):
    s = se_r2(r2, n)
    return max(0.0, r2 - 1.96 * s), r2 + 1.96 * s

trainN = read_training_N()
rows = []
nimp = 0
for idd, name in TRAITS:
    g, l, n = read_r2(idd)
    delta = l - g
    pct = 100 * delta / g if g else float("nan")
    imp = l > g
    nimp += int(imp)
    glo, ghi = ci(g, n); llo, lhi = ci(l, n)
    rows.append([name, idd,
                 f"{g:.4f}", f"{glo:.4f}-{ghi:.4f}",
                 f"{l:.4f}", f"{llo:.4f}-{lhi:.4f}",
                 f"{delta:+.4f}", f"{pct:+.1f}",
                 "yes" if imp else "no", int(n),
                 trainN.get(idd, "")])

hdr = ["Trait", "AoU phenotype ID",
       "EUR R2 (Global / All-by-All)", "Global 95% CI",
       "EUR R2 (Local / FELIXassoc EUR-specific)", "Local 95% CI",
       "delta R2 (Local - Global)", "relative change (%)",
       "Local improves", "UKB EUR validation N", "AoU EUR training N (local)"]

out_tsv = os.path.join(PRS_DIR, "SupplementaryTable9_PRS.tsv")
with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(rows)
print("wrote", out_tsv)
print(f"Local > Global (EUR) in {nimp}/{len(TRAITS)} traits")

# optional xlsx
try:
    import openpyxl
    from openpyxl.styles import Font
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Stable9_PRS"
    ws.append(["Supplementary Table 9. European polygenic score prediction R2: "
               "FELIXassoc EUR-specific (Local) vs All-by-All EUR (Global) summary "
               "statistics, UK Biobank British-ancestry validation. Both scores are "
               "evaluated in the same individuals, so the Local-vs-Global comparison is "
               "internally controlled. 95% CI from the large-sample SE 2*sqrt(R2)*(1-R2)/sqrt(N)."])
    ws.append([]); ws.append(hdr)
    for r in rows: ws.append(r)
    for c in ws[3]: c.font = Font(bold=True)
    wb.save(os.path.join(PRS_DIR, "SupplementaryTable9_PRS.xlsx"))
    print("wrote", os.path.join(PRS_DIR, "SupplementaryTable9_PRS.xlsx"))
except Exception as e:
    print("(xlsx skipped:", e, ")")
