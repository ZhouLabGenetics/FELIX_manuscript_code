# FELIX manuscript — analysis code

Code accompanying the FELIX manuscript. FELIX (Full-cohort Efficient Local
ancestry-Integrated miXed-model framework) has two components: **FELIXla**, a
haplotype-resolved storage format, and **FELIXassoc**, the association module.
The method is implemented in the FELIX R package:
https://github.com/LinfengHu-1/SAIGE-Tractor

## Layout

- `aou/` — commands run on the All of Us Researcher Workbench: the FELIXassoc
  null-model fit (step 1), the FELIXla conversion, and the score tests (step 2).
  The real-data analysis is a direct application of the FELIX package; these are
  the exact invocations.
- `simulation/` — three-way admixture simulation for type I error, power and the
  Tractor-Mix / SAIGE benchmarks. `R/` holds the generators and analysis,
  `slurm/` the job scripts (numbered in run order).
- `prs/` — polygenic score pipeline: prep sumstats → PRS-CS → intersect → PLINK2
  scoring → R² evaluation in UK Biobank, plus the Figure 4 and Supplementary
  Table 9 scripts. `run_phenotype.sh` chains the per-phenotype steps.
- `figures/` — scripts for the manuscript figures 

Absolute paths reflect the cluster and workbench the scripts were run on.
