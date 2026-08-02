# Runbook — UKB multi-ancestry PRS (new cluster / UGER / qsub)

Living document. Records the goal, the design, every input, and the exact commands so the whole analysis can be replicated. Update it whenever the pipeline changes.

## 0. Goal & design

**Question (Ying + Wei):** given a full-cohort FELIXassoc GWAS, does its full-cohort effect estimate (`BETA_c_ancALL`) predict UKB phenotypes better than conventional GWAS summary statistics? We test this with **P+T first (fast go/no-go)**, then PRS-CS only if P+T is promising (esp. the modest EUR gain they expect).

**3 predictors** (score = SNP A1 A2 BETA P; SNP = `CHR:POS`, A1 = effect allele = Allele2): \| name \| source \| type \| LD panel for clumping \| \|---\|---\|---\|---\| \| `meta` \| ABA `<id>_META.tsv` `BETA` \| multi-ancestry \| eur10k **and** prop10k \| \| `felix` \| FELIX `merged_FELIX_*` `BETA_c_ancALL` \| multi-ancestry \| eur10k **and** prop10k \| \| `matched_<ANC>` \| ABA `<id>_<ANC>.tsv` (EUR/AFR/EAS/SAS→CSA) \| single-ancestry \| ancestry-matched \|

**4 validation sets (UKB):** EUR, AFR, EAS, CSA. **11 quantitative traits.** **Build:** predictors (ABA/FELIX from AoU) are **hg38** (All of Us is GRCh38); the UKB bgen here AND the PRS-CS snpinfo/LD panels are **hg19**. Confirmed by rs11240767 = 728951 in the bgen (hg19) vs 793571 in the predictors (hg38). So predictors are lifted **hg38→hg19** once (`liftover_predictors.sh`), and genotype variant IDs are forced to `CHR:POS` (`plink2 --set-all-var-ids @:#`) — everything ends up hg19 CHR:POS and matches.

**Ancestry-matched coverage (predictor 2 only; 1 and 3 exist for all 11 × 4):** EUR/AFR all 11; EAS 6 (3007070, 3009744, 3022192, 3024929, 3027114, height); CSA (ABA SAS) 8 (3006923, 3007070, 3009744, 3013721, 3022192, 3024929, 3035995, BMI).

**LD panels (for clumping only):** `eur10k` = 10k random UKB-EUR; `prop10k` = 10k mixed, proportional to the pooled-discovery ancestry sizes from height's FELIX `N_haplo` (anc1 AFR 98909, anc2 EAS 13711, anc3 EUR 310264, anc4 AMR 33639, anc5 SAS 13319); AMR dropped (not in UKB), renormalized over EUR/AFR/EAS/CSA → `EUR:7113 AFR:2268 EAS:314 CSA:305`. Single-ancestry predictor uses its matched panel. eur10k-vs-prop10k = the "is EUR LD good enough" sensitivity.

## 1. Cluster facts

-   Scheduler UGER/`qsub`; only `apptainer` + manually-installed tools in `/humgen/atgu1/fin/lhu/tools`.
-   `plink2` = `tools/plink2` (v2.00a6LM) — used for bgen→bed extraction.
-   `plink` = `tools/plink` (v1.9.0-b.7.11) — used for `--clump` and `--score` (validated dialect).
-   R = Broad dotkit **R-4.1** (`set +u; source /broad/software/scripts/useuse; use R-4.1; set -u`), base R only.
-   Genotypes: `/broad/ukbb/imputed_v3/ukb_imp_chr${c}_v3.bgen` (+ `.bgi`), **hg19** (rs11240767=728951).
-   bgen sample file: `/humgen/atgu1/fin/wzhou/projects/survival_analysis/realdata/UKBB/geno/ukb31063.autosomes.sample.rmFirstTwoLines`.
-   Per-ancestry sample lists: `/humgen/atgu1/fin/wzhou/projects/UKB/genotype/ukb.<ANC>.for_grm.*.pruned.plink.fam` (EUR fam has `.1e7.`).
-   HM3 list: `tools/PRS_ref/ldblk_ukbb_eur/snpinfo_ukbb_hm3`.
-   Master dir: `/humgen/atgu1/fin/lhu/projects/saige_tractor/prs`.

## 2. What to transfer (nothing huge)

Build the score files LOCALLY, upload only the small outputs.

``` bash
# LOCAL (mac):
cd ~/Desktop/5pop_st/prs && bash build_inputs_local.sh          # -> predictor_inputs/
# upload:
DST=lhu@<newcluster>:/humgen/atgu1/fin/lhu/projects/saige_tractor/prs
rsync -avz ~/Desktop/5pop_st/prs/predictor_inputs/            ${DST}/predictor_inputs/
rsync -avz ~/Desktop/5pop_st/prs/newcluster/                  ${DST}/scripts/
rsync -avz ~/Desktop/5pop_st/prs/bootstrap_r2.R               ${DST}/scripts/
rsync -avz ~/Desktop/5pop_st/prs/subset_ukbb_biomarkers.tsv   ${DST}/pheno/
rsync -avz ~/Desktop/5pop_st/prs/final_samples.txt.bgz        ${DST}/covar/
rsync -avz ~/Desktop/5pop_st/prs/sample_id_mapping.txt        ${DST}/
```

Already on the cluster (no transfer): genotypes, bgen sample file, `PRS_ref/` LD panels, plink2, plink1.9, saigeqtl sif.

## 3. Pipeline (run on the cluster, from the scripts dir)

``` bash
cd /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/scripts
# config.sh is already filled (BGEN_SAMPLE, PROP, plink paths). Review once.

rm -rf /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/geno_bed/*   # clear the failed run's temps
bash 00_setup_dirs.sh                       # make the tree (now: geno_bed/{full,ldset,eurval})
bash 01_prep_refs.sh                        # HM3 lists + keeps incl. ldset/eurval/needed

# genotype extraction: 22 tasks (one per chr). bgenix pulls HM3-only variants from the
# bgen via its .bgi index (fast; confirms the rsID list works) -> small HM3 bgen; then
# plink2 splits it into ldset + eurval (FULL EUR) beds, renaming to CHR:POS. Memory capped.
SGE_TASK_ID=22 bash 02_geno_prep.qsub       # SMOKE TEST (chr22): watch the HM3 variant count (should be ~tens of k)
qsub 02_geno_prep.qsub                      # 22 tasks in parallel

# SANITY after chr22: predictor CHR:POS should now overlap the genotypes heavily (hg19)
awk 'NR==FNR{g[$2];next} FNR>1{t++; if($1 in g)o++} END{printf "overlap %d/%d = %.1f%%\n",o,t,100*o/t}' \
  geno_bed/ldset/chr22.bim predictor_inputs/BMI/meta.txt   # expect a high %, not 0.7%

# BUILD FIX: AoU predictors are hg38, UKB bgen + PRS-CS refs are hg19. Lift predictors
# hg38->hg19 (rs11240767: bgen 728951 hg19 vs predictor 793571 hg38). Run once.
bash liftover_predictors.sh                 # backs up to predictor_inputs_hg38_bak/, rewrites keys to hg19

# P+T: one array task per trait
SGE_TASK_ID=6 bash 03_pt.qsub               # SMOKE TEST trait #6 = 3024929 (platelets)
qsub 03_pt.qsub                             # 5e-8
qsub -v CLUMP_P1=5e-6 03_pt.qsub            # 5e-6 sensitivity

# incremental R2 + bootstrap — PARALLEL by trait, PHASED so partial results survive walltime.
# Each task: (1) EUR-10k -> (2) AFR/EAS/CSA -> (3) full EUR (last). R = R-4.1 via useuse.
qsub 04_bootstrap.qsub                       # 5e-8, B=2000 (phased)
qsub -v CLUMP_P1=5e-6 04_bootstrap.qsub      # 5e-6
qsub -v B=1000 04_bootstrap.qsub             # faster full-EUR tail (looser CIs)
# main EUR incR2_<thr>.tsv = 10k fallback until full EUR finishes (tell them apart by the N
# column); labeled 10k copy always kept at incR2_<thr>_eur10k.tsv. 04b is now redundant.

# aggregate + plot for the boss (base R; no ggplot dep)
use Python-3.9
python3 aggregate_results.py  /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/results  5e-8
set +eu; source /broad/software/scripts/useuse; use R-4.1; set -eu    # nonzero when already loaded
Rscript plot_results.R  /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/results  5e-8

# PRS-CS PREP — safe to run NOW in parallel with the bgen extraction
bash prscs_prep.sh
```

## 4. Predictors produced per (trait × validation ancestry)

Built by `make_predictor_inputs.py` (score files) → `pt_predict.sh`/`run_pt.sh` (PRS). \| score name \| meaning \| LD panel \| role \| \|---\|---\|---\|---\| \| `felix__eur10k`, `felix__prop10k` \| FELIX full-cohort `BETA_c_ancALL` \| eur10k / prop10k \| **main test** \| \| `meta__eur10k`, `meta__prop10k` \| AllxAll multi-ancestry META \| eur10k / prop10k \| main benchmark \| \| `matched_<VAL>__<VAL>` \| AllxAll ancestry-matched (global) \| matched \| main benchmark \| \| `felix_tract_<VAL>__<VAL>` \| FELIX ancestry-TRACT `BETA_c_anc<k>` (local) \| matched \| supp. global-vs-local \|

`k`: AFR=1, EAS=2, EUR=3, CSA=5 (SAS). So per cell you get up to **6 PRS**.

## 5. Outputs & how to read them

-   `geno_bed/ldset/chr<c>.{bed,bim,fam}` — HM3, CHR:POS ids, all LD panels + AFR/EAS/CSA targets (\~35k people; subset at run time via `--keep keeps/<panel>.keep`). `geno_bed/eurval/chr<c>` — the EUR validation target (FULL 451,509; `EUR_VALID_N=""`).
-   `results/<trait>/<VAL>/<score>_<thresh>.sscore` — per-PRS scores.
-   `results/<trait>/<VAL>/incR2_<thresh>.tsv` (+ `_diff.tsv`) — incremental R² + 95% CI per PRS; `_diff` = each PRS **minus FELIX** (baseline `felix__eur10k`): negative with CI below 0 ⇒ FELIX wins; also holds `felix__prop10k − felix__eur10k` = LD sensitivity.
-   `results/summary_<thresh>.tsv`, `results/diff_vs_felix_<thresh>.tsv` — tidy master tables (all trait×ancestry×PRS) from `aggregate_results.py`.
-   `results/fig_main_predictors_<thresh>.pdf` — FELIX vs META vs matched, faceted by ancestry.
-   `results/fig_ld_sensitivity_<thresh>.pdf` — eur10k vs prop10k (points near the y=x line ⇒ LD choice small).
-   `results/fig_tract_global_local_<thresh>.pdf` — global (AllxAll) vs local (FELIX tract) per ancestry.

**Two claims for the greenlight:** (a) `felix` ≥ `meta` and `matched` across traits/ancestries (esp. modest EUR gain) — read `fig_main_predictors` + `diff_vs_felix`; (b) eur10k ≈ prop10k for the multi-ancestry predictors — read `fig_ld_sensitivity`.

## 5b. EUR global-vs-local (tract) — the old AFR analysis, now for EUR (+ all 4)

Same design as the old-cluster AFR global-vs-local, on the new UKB cohorts. **global** = AllxAll ancestry-matched (`matched_<VAL>`), **local** = FELIX ancestry-tract (`felix_tract_<VAL>`, EUR = `BETA_c_anc3`). Both are produced automatically by `run_pt.sh` for all 4 ancestries and land in the same `incR2_<thresh>.tsv`; the comparison plot is `fig_tract_global_local_<thresh>.pdf`. EUR is the panel of interest (11 traits), the others are supplementary. No extra run needed — it's folded into the main P+T.

## 6. Next stage (only if P+T is promising): PRS-CS

-   **Prep (`prscs_prep.sh`):** writes rsID-keyed sumstats (`prscs_inputs/<trait>/<pred>.txt`, `SNP A1 A2 BETA P`) — PRS-CS matches by rsID. Maps predictor CHR:POS -> rsID via `snpinfo_ukbb_hm3`. **Run only after the predictors are lifted to hg19** (snpinfo is hg19, so hg19 predictors map; hg38 would not). No build-specific code change needed — the lift handles it. The script prints `SANITY BMI/meta: mapped X/Y` — expect X ≈ Y.
-   LD panels already at `tools/PRS_ref/ldblk_ukbb_{eur,afr,eas,sas}` (HM3). Panel↔validation: eur→EUR, afr→AFR, eas→EAS, sas→CSA. Multi-ancestry predictors (felix, meta): use the EUR panel first, then the proportional-mix sensitivity.
-   **Still needed before running PRS-CS:** (1) PRScs software (`git clone` PRScs; needs python with `scipy` + `h5py` — check if the saige sif has them, else a small env);
    (2) per-trait `--n_gwas` (approx from FELIX `N_haplo_ancALL`/2). Scores + bootstrap reuse the same target beds + `run_bootstrap.sh`.

## Change log

-   2026-07-30: initial pipeline (bgen→HM3 beds, plink1.9 P+T on eur10k/prop10k/matched, bootstrap incR² vs FELIX). Genotypes switched from the pruned for_grm files to HM3-from-imputed-bgen per Wei/Ying. Fixed 02 output-dir-missing bug.
-   2026-07-30: two smoke-test fixes. (1) renamed the `GROUPS` array → `GRP` in `02_geno_prep.qsub` — `GROUPS` is a special bash variable (user group IDs), so assignments were ignored and `${GROUPS[0]}` returned the GID `1015`. (2) plink2 `--sample` needs the 2-line Oxford header; the `...rmFirstTwoLines` file has it stripped, so `01_prep_refs.sh` now rebuilds `SAMPLE_PLINK=${BASE}/ukb_imp.sample` (ID_1 ID_2 missing + type row + the data rows) and `02` uses it. **Verify** the smoke test now prints `group=EUR` and keeps a non-zero sample count (\~451509 for EUR); if 0 kept, the for_grm `.fam` FID/IID don't match the `.sample` ID_1/ID_2 — tell me.
-   2026-07-30: the `...rmFirstTwoLines` file turned out to be a SINGLE column of sample IDs, so `01_prep_refs.sh` now writes `ID_1 ID_2 missing` = `id id 0` per sample.
-   2026-07-30: qsub copies the job script to a spool dir, so `$(dirname "$0")/config.sh` wasn't found (BASE empty -\> `/geno_bed/...` permission errors). All scripts now `source` config by ABSOLUTE path and call siblings via `${SCRIPTS_DIR}` (=`${BASE}/scripts`). **The scripts MUST live at `${BASE}/scripts`** for this to work.
-   2026-07-30: added the FELIX ancestry-tract predictor (`felix_tract_<ANC>` = `BETA_c_anc<k>`) to `make_predictor_inputs.py` + `run_pt.sh` → EUR global-vs-local (and all 4) folded into the main P+T (§5b). Added `aggregate_results.py`, `plot_results.R` (3 figures for the boss), and `prscs_prep.sh` (rsID-keyed PRS-CS inputs, safe to run now).
-   2026-07-30: geno-prep OOM + disk blowup fixed. Cause: 132 tasks each re-imported the bgen and `--make-bed` spilled a full-cohort temporary pgen (\~800 GB total); plink2 also reserved half the node RAM and got cgroup-OOM-killed. New design: 22 tasks (one per chr), `--make-pgen --keep needed` (filters samples on write, no full temp) → cheap subsets to `ldset`/`eurval` → delete pgen; `--memory ${PLINK2_MEM}`; submit with `-tc 4`. Clumping panels are all ≤10k (EUR-matched uses `eur10k`); only EUR *scoring* needs many people, and that target is subsampled to `EUR_VALID_N=50000` for the sanity check (set "" for full). pt_predict/run_pt now read `ldset`(+`--keep`)/`eurval` instead of per-group beds.
-   2026-07-30: BUILD-MISMATCH fix. `hm3/chr*.extract` (from snpinfo BP) overlapped the hg19 predictors only 0.7% → snpinfo BP is NOT hg19. Now select HM3 by **rsID** (`hm3/hm3_rsids.txt`, snpinfo `SNP` col) at import, then `--set-all-var-ids @:#` in the subset step to get CHR:POS from the bgen's own hg19 positions. No liftOver/downloads needed. Sanity-check the overlap (command above) after chr22 — expect a high %. If the import keeps almost NO variants, the bgen variant IDs aren't rsIDs and we'd fall back to position-based selection.
-   2026-07-30: full EUR (451,509) restored (`EUR_VALID_N=""`). To keep it affordable, `02` now uses **bgenix** to pre-extract HM3-only variants from the bgen (index-based, fast, and confirms the rsID list) into a small HM3 bgen, then plink2 splits that into ldset/eurval beds (`--set-all-var-ids @:#`). Intermediate HM3 bgen is deleted per chr. `prop10k` confirmed = AoU-discovery proportions (EUR 71/AFR 23/EAS 3/CSA 3, ex-AMR from height N_haplo), i.e. mimic AoU, NOT UKB.
-   2026-07-30: `--keep` matched 0 samples. Our `.sample` has ID_1=ID_2=eid (FID=eid), but the keeps were `$1,$2` from the for_grm .fam = FID=0. plink2 matches FID AND IID → no match. Fixed: `01` builds keeps as `$2,$2` (FID=IID=eid). bgenix HM3 extraction itself works (16,360 variants on chr22). Re-run `01` to rebuild all keep files, then re-smoke `02`.
-   2026-07-30: **BUILD MISMATCH found.** Genotype bim IDs are correct CHR:POS, but chr22 position overlap with hg19 predictors was 11/16360 — below random =\> different builds. Predictors are hg19 (META 1:793571 = rs11240767 GRCh37); the /broad/ukbb/imputed_v3 bgen (and PRS-CS snpinfo/LD panels) are **hg38**. Confirm bgen build: `bgenix -g ukb_imp_chr1_v3.bgen -incl-rsids <(echo rs11240767) -list` (858191=hg38, 793571=hg19). Fix = `liftover_predictors.sh` lifts predictor CHR:POS hg19-\>hg38 (UCSC liftOver + hg19ToHg38 chain, auto-downloaded; backs up to predictor_inputs_hg19_bak/). This ALSO makes prscs_prep.sh correct (snpinfo is hg38). Run once before `03`.
-   2026-07-30: **direction corrected.** bgenix -list gave rs11240767 = 728951 in the UKB bgen (=hg19), vs 793571 in the predictors (=hg38, since AoU is GRCh38). So UKB bgen + PRS-CS snpinfo/LD panels are hg19, predictors are hg38. liftover_predictors.sh lifts predictors **hg38-\>hg19** (hg38ToHg19 chain), matching genotypes AND the hg19 PRS-CS refs. (My earlier "predictors hg19 / bgen hg38" was backwards — the bgen is authoritative for hg19.)
-   2026-07-30: liftover done LOCALLY (cluster installs are hard). `liftover_predictors_local.py` (pyliftover + hg38ToHg19 chain in \~/Desktop/5pop_st/prs/lift/) lifted predictor_inputs hg38-\>hg19 in place (backup predictor_inputs_hg38_bak/); validated 1:793571-\>1:728951 (matches bgen). 100% kept (402/1.04M unmapped). `build_felix_tract_hg19.py` builds the felix_tract\_<ANC> files directly in hg19 (were missing from the earlier build). After both, re-upload predictor_inputs/ to the cluster. Do NOT re-run build_inputs_local.sh (it writes hg38).
- 2026-07-30: **snpinfo build clarified = hg19** (earlier change-log lines calling it hg38
  are superseded). With predictors lifted to hg19, `prscs_prep.sh` needs NO change: it maps
  predictor CHR:POS(hg19) -> rsID via snpinfo(hg19). Added a `SANITY BMI/meta: mapped X/Y`
  print — X≈Y confirms alignment; near-0 would mean predictors weren't lifted or snpinfo
  isn't hg19. Run prscs_prep AFTER re-uploading the hg19 predictor_inputs.
- 2026-07-30: bootstrap parallelized. EUR validation is ~400k people, so B*fit is slow
  (>15 min/trait serial). `run_bootstrap.sh` now takes an optional single-trait arg and
  `04_bootstrap.qsub` runs one trait per array task (1-11) -> ~one trait's wall time.
  `B` is an env knob (default 2000; use 1000 to halve). R switched from the saigeqtl sif
  to Broad dotkit **R-4.1** (`set +u; source /broad/software/scripts/useuse; use R-4.1; set -u`);
  config RUN_R="Rscript". covars.tsv + geno2eid map are written atomically (temp+mv) so
  parallel tasks don't corrupt them. plot_results.R also runs under R-4.1.
- 2026-07-30: added EUR-10k INSURANCE bootstrap (boss OK'd a 10k EUR subset, wider CI, keep
  B=2000). run_bootstrap.sh gained env knobs EUR_KEEP (--keep EUR to a subset), ONLY_VAL
  (restrict ancestries), OUTSUF (output suffix). `04b_bootstrap_eur10k.qsub` runs EUR-only
  on keeps/eur10kval.keep (10k EUR, disjoint from the eur10k LD panel, seed 45) -> writes
  results/<trait>/EUR/incR2_<thr>_eur10k.tsv WITHOUT touching the full-EUR incR2_<thr>.tsv.
  Uses the existing full-EUR sscores (just --keep'd), so it's seconds/trait. Safe to run
  alongside 04_bootstrap.qsub. If full EUR times out, use the _eur10k EUR rows instead.
- 2026-07-30: `04_bootstrap.qsub` rewritten PHASED (per trait): phase 1 EUR-10k (writes
  incR2_<thr>_eur10k.tsv AND copies it to the main EUR incR2_<thr>.tsv as a fallback),
  phase 2 AFR/EAS/CSA, phase 3 FULL EUR last (overwrites the main EUR file if it finishes).
  So a walltime kill in phase 3 still leaves all 4 ancestries + FELIX comparison (EUR@10k).
  Distinguish 10k vs full EUR by the N column. Uses run_bootstrap.sh's ONLY_VAL/EUR_KEEP/
  OUTSUF knobs. 04b_bootstrap_eur10k.qsub is now redundant (kept for standalone use).
  aggregate_results.py picks up whatever is in the main incR2_<thr>.tsv (10k or full).

## 7. PRS-CS pipeline (main manuscript result)

Design: per validation ancestry T, run PRS-CS for each predictor (felix ancALL, meta, matched_T)
re-weighted by T's UKBB LD panel (EUR->eur, AFR->afr, EAS->eas, CSA->sas), score on the T
validation set, bootstrap incremental R2 (baseline = felix). Same 4 targets, 11 traits.

Uses: python3.9 (scipy,h5py) + R-4.1 via `use`; PRScs at `tools/PRScs/PRScs.py`;
LD panels `tools/PRS_ref/ldblk_ukbb_{eur,afr,eas,sas}`; n_gwas from `n_gwas_prscs.tsv`
(felix=TOTAL_local, meta=TOTAL_global, matched_<T>=<T>_global; built locally, do NOT parse
the ragged sample_size_table on the cluster).

### Upload (in addition to §2)
```
DST=lhu@<cluster>:/humgen/atgu1/fin/lhu/projects/saige_tractor/prs
rsync -avz ~/Desktop/5pop_st/prs/n_gwas_prscs.tsv   ${DST}/           # -> ${BASE}/n_gwas_prscs.tsv
rsync -avz ~/Desktop/5pop_st/prs/newcluster/        ${DST}/scripts/   # updated scripts
# (PRScs already git-cloned at tools/PRScs; LD panels already in tools/PRS_ref)
```

### Run
```bash
cd /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/scripts
bash prscs_prep.sh                                  # rsID sumstats + prscs_hm3.bim + prscs_combos.txt
NC=$(wc -l < /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/prscs_combos.txt)

# PRS-CS is slow (~3-4 h for the biggest chromosomes at 1000 MCMC iters) -> parallelize PER
# CHROMOSOME (combo x chr). h_rt=12h/task (max chr ~4h; drop to 8h if your queue caps lower).
B=/humgen/atgu1/fin/lhu/projects/saige_tractor/prs
SGE_TASK_ID=1 bash 05_prscs.qsub                    # SMOKE TEST one chr; confirm posterior name:
ls ${B}/prscs_out/*/*/*_pst_eff_*_chr*.txt | head

# --- Alkaline Phosphatase (3035995) FIRST (boss wants the PRS-CS trend validated on it) ---
grep '^3035995 ' ${B}/prscs_combos.txt > ${B}/prscs_combos_alp.txt
NA=$(wc -l < ${B}/prscs_combos_alp.txt)             # ~11 combos -> 242 chr-tasks
qsub -t 1-$((NA*22)) -tc 300 -v COMBOS=${B}/prscs_combos_alp.txt 05_prscs.qsub
qsub -t 1-${NA}      -v COMBOS=${B}/prscs_combos_alp.txt 06_prscs_score.qsub   # after 05 (ALP) done
qsub -t 9-9          07_prscs_bootstrap.qsub        # ALP is trait #9 (DUMPs draws by default)

# --- the other 10 traits ---
grep -v '^3035995 ' ${B}/prscs_combos.txt > ${B}/prscs_combos_rest.txt
NR=$(wc -l < ${B}/prscs_combos_rest.txt)
qsub -t 1-$((NR*22)) -tc 300 -v COMBOS=${B}/prscs_combos_rest.txt 05_prscs.qsub
qsub -t 1-${NR}      -v COMBOS=${B}/prscs_combos_rest.txt 06_prscs_score.qsub
qsub                 07_prscs_bootstrap.qsub        # all 11 (re-does ALP harmlessly)

# aggregate + main figure (base R / python3.9)
use Python-3.9
python3 aggregate_results.py /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/results PRSCS
# -> summary_PRSCS.tsv, diff_vs_felix_PRSCS.tsv  (main figure/table scripts: pending violin-layout decision)
```

### Outputs
- `prscs_out/<trait>/<T>/<pred>_pst_eff_..._chr<c>.txt` — per-chr posterior effects (rsID).
- `results/<trait>/<T>/<pred>_PRSCS.sscore` — PRS-CS scores on the target.
- `results/<trait>/<T>/incR2_PRSCS.tsv` (+ `_diff`, `_eur10k`) — incR2 + bootstrap, baseline felix.
- `results/summary_PRSCS.tsv`, `results/diff_vs_felix_PRSCS.tsv` — master tables for the manuscript.

### Notes
- 05 is ONE task per combo (~124 tasks total), each running all 22 chr internally (not NC*22).
  Out-of-range tasks skip. Optionally cap concurrency with `-tc` (e.g. `qsub -t 1-${NC} -tc 200`).
  Each task ~2-6 h (22 chr sequential); if that risks walltime, split with CHR (per-chr array).
- If PRScs writes a different posterior filename (version/phi), fix the pattern in 05's skip
  check and 06's `POST=` glob to match `ls prscs_out/*/*/*_pst_eff_*`.
- Figures/colors: FELIX ancALL=#AA3377, All by All META=#0077BB, All by All ancestrally-matched
  =khaki2 (#EEE685); short-dash legends; hi-res PDF. (Same spec as the P+T violin backup.)

## Change log (cont.)
- 2026-07-31: legend renamed -> "Global (All by All)" / "Local (FELIXassoc)" in both
  plot_partial_prs.R (AFR) and plot_eur_global_local.R (EUR) + standalone legend PNG.
  Built PRS-CS pipeline: n_gwas_prscs.tsv (local, from sample_size_table), config PRS-CS
  block, prscs_prep.sh (+prscs_hm3.bim +prscs_combos.txt), 05_prscs.qsub (posteriors,
  combo x chr), 06_prscs_score.qsub (score+combine), 07_prscs_bootstrap.qsub (phased,
  BASELINE=felix, tag PRSCS). run_bootstrap.sh gained a BASELINE env. R=R-4.1, py=Python-3.9.

## 8. Figures (P+T backup + PRS-CS main)

Colors (from scripts_R/fig05_extra.R): FELIX ancALL=#AA3377, All by All META=#0077BB,
All by All ancestrally-matched=khaki2 #EEE685. Short-dash labels. Hi-res PDF+PNG.

Two violin layouts (scripts take <results_dir> <tag> [outstem]; tag=5e-8 for P+T, PRSCS for PRS-CS):
- **Layout A (by method)** `plot_violin_bymethod.R` — x=3 methods, violin over 11 traits + trait
  dots, 4 ancestry panels, shared y. Reads `summary_<tag>.tsv` (NO re-run).
- **Layout B (by trait)** `plot_violin_bytrait.R` — x=11 trait names, 3 dodged violins of the
  2000 bootstrap draws, 4 panels, shared y. Needs `*_draws.tsv.gz` from a bootstrap run with DUMP=1.

Generate (base R via R-4.1):
```bash
set +eu; source /broad/software/scripts/useuse; use R-4.1; set -eu
cd /humgen/atgu1/fin/lhu/projects/saige_tractor/prs/results
# P+T (already have summary_5e-8): layout A now; layout B needs draws (see below)
Rscript ../scripts/plot_violin_bymethod.R . 5e-8  fig_violin_bymethod_5e-8
Rscript ../scripts/plot_violin_bytrait.R  . 5e-8  fig_violin_bytrait_5e-8
# PRS-CS (after 07): both
Rscript ../scripts/plot_violin_bymethod.R . PRSCS fig_violin_bymethod_PRSCS
Rscript ../scripts/plot_violin_bytrait.R  . PRSCS fig_violin_bytrait_PRSCS
```

**To get layout B (bootstrap-draw violins) you must bootstrap with DUMP=1:**
```bash
qsub -v DUMP=1 04_bootstrap.qsub                 # P+T: re-run to write *_draws.tsv.gz
qsub 07_prscs_bootstrap.qsub                     # PRS-CS: DUMPs draws BY DEFAULT (no -v needed)
```
(DUMP writes `results/<trait>/<VAL>/incR2_<tag>_draws.tsv.gz`; layout B falls back to the
eur10k EUR draws if full-EUR draws are absent.)

Global-vs-local supplementary (2-method, blue/pink): `plot_partial_prs.R` (AFR) and
`plot_eur_global_local.R` (EUR) + standalone `legend_global_local.png`.

## Change log (cont.)
- 2026-07-31: bootstrap_r2.R gained `--dump_boot` (writes B x score draws); run_bootstrap.sh
  gained DUMP env. Added plot_violin_bymethod.R (layout A, from summary, no re-run) and
  plot_violin_bytrait.R (layout B, from *_draws.tsv.gz, needs DUMP=1). Both for P+T + PRS-CS,
  fig05 colors, shared y, short-dash labels. P+T layout A generated locally.
- 2026-08-01: FIX exit_status 141 (SIGPIPE). run_bootstrap.sh picked a representative sscore
  with `ls .../*_TAG.sscore | head -1`; `head` closes the pipe -> `ls` SIGPIPE=141, and under
  set -e/pipefail the script died silently (empty .e, .o cut at the PHASE-1 echo, task gone).
  Racy (buffer timing) so only some array tasks failed; worsened by the generic `*_TAG.sscore`
  glob matching many files. Fixed pipe-free (a for-loop that breaks on the first match). Not a
  memory/walltime/data issue (sscores were complete at 451,510 lines).
