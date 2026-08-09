# FELIX manuscript analysis scripts

Analysis code accompanying the FELIX manuscript. FELIX (Full-cohort Efficient Local
ancestry-Integrated miXed-model framework) has two components: **FELIXla**, a
haplotype-resolved storage format, and **FELIXassoc**, the association module.
The method is implemented in the FELIX R package:
[https://github.com/ZhouLabGenetics/FELIX](https://github.com/ZhouLabGenetics/FELIX.git)

## Layout

- `aou/` — example commands run on the All of Us Researcher Workbench: the FELIXla conversion,
  the FELIXassoc null-model fit (step 1), and the score tests (step 2).
  The real-data analysis is a direct application of the FELIX package.
- `simulation/` — three-way admixture simulation for FELIX type I error, power and 
  benchmarking. `R/` holds the generators and analysis,
  `slurm/` the job scripts.
- `prs/` — polygenic score pipeline: prep sumstats → PRS-CS → intersect → PLINK2
  scoring → R² evaluation in UK Biobank, plus the Figure 4 and Supplementary
  Table 9 scripts. `run_phenotype.sh` chains the per-phenotype steps.
- `figures/` — scripts for the manuscript figures 

Generated results and controlled-access inputs are not included.
