# Multi-ancestry PRS pipeline (P+T and PRS-CS)

Polygenic score evaluation used in the FELIX manuscript. Discovery summary statistics are from
All of Us (FELIXassoc full-cohort `BETA_c_ancALL`; All by All multi-ancestry meta; All by All
ancestry-matched). Validation is in four Pan-UK Biobank cohorts (EUR, AFR, EAS, CSA) across 11
quantitative traits. Prediction is quantified as incremental R2 with a 2,000-replicate bootstrap.
Paths in the scripts are for the Broad UGER cluster; adapt `config.sh` for other systems. Full
step-by-step notes are in `RUNBOOK_newcluster_prs.md`.

## Layout / order
- `config.sh` shared paths and helpers (edit once).
- `make_predictor_inputs.py`, `build_inputs_local.sh` build score files (SNP A1 A2 BETA P) from the
  AoU summary statistics (run locally; predictors are hg38 and are lifted to hg19).
- `liftover_predictors.sh` / `liftover_predictors_local.py` lift predictor keys hg38 to hg19.
- `00_setup_dirs.sh`, `01_prep_refs.sh` directory tree, HapMap3 rsID list, per-ancestry keep lists,
  LD panels (`eur10k`, `prop10k`).
- `02_geno_prep.qsub` extract HapMap3 variants from UK Biobank imputed bgen (bgenix + plink2).
- P+T: `pt_predict.sh`, `run_pt.sh`, `03_pt.qsub` (clump + score); `04_bootstrap.qsub` /
  `run_bootstrap.sh` (incremental R2 + bootstrap, phased so EUR-10k precedes full EUR).
- PRS-CS: `prscs_prep.sh` (rsID sumstats + bim + combos), `05_prscs.qsub` (posteriors, per combo x
  chromosome), `06_prscs_score.qsub` (score + combine), `07_prscs_bootstrap.qsub` (incremental R2 +
  bootstrap, baseline FELIX). `n_gwas_prscs.tsv` holds the per-trait GWAS sample sizes.
- Aggregate / tables: `aggregate_results.py`, `make_prscs_supp_tables.py`,
  `make_prscs_combined_table.py`.
- Figures: `plot_violin_bytrait.R`, `plot_violin_bymethod.R` (main; FELIX vs META vs
  ancestry-matched), `plot_global_local_supp.R`, `plot_eur_global_local.R` (supplementary
  local vs global), `plot_results.R`.
- `bootstrap_r2.R` incremental R2 with paired bootstrap CIs (`--dump_boot` writes per-replicate
  draws for the by-trait violins).

## Reproduce (PRS-CS)
```
bash prscs_prep.sh
NC=$(wc -l < prscs_combos.txt)
qsub -t 1-$((NC*22)) -tc 300 05_prscs.qsub     # posteriors (combo x chromosome)
qsub -t 1-${NC}      06_prscs_score.qsub        # score + combine
qsub                 07_prscs_bootstrap.qsub    # incremental R2 + bootstrap (dumps draws)
python3 aggregate_results.py results PRSCS
Rscript plot_violin_bytrait.R  results PRSCS out/fig4a
Rscript plot_violin_bymethod.R results PRSCS out/fig4b
python3 make_prscs_combined_table.py results PRSCS
```
