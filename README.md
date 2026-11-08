# Stair-climbing frequency and COPD: two-sample Mendelian randomization

Analysis code and frozen summary-level inputs for:

> Lou J, Yu L, Xiao L, Shi Z. Associations of genetically proxied stair-climbing frequency with chronic obstructive pulmonary disease-related outcomes and lung function: a two-sample Mendelian randomization study. *Manuscript submitted for publication.*

The scripts reproduce every univariable MR estimate, sensitivity analysis, supplementary analysis, table source file and Figure 1 reported in the manuscript and its Supporting Information. All inputs are publicly available GWAS summary statistics; no individual-level data are included.

## Repository layout

| Folder | Contents |
|---|---|
| `inputs/` | Frozen, harmonized summary-level inputs read by the scripts (exposure–outcome SNP associations, negative-control and MVMR extracts, GWAS metadata, instrument annotation, dataset descriptors) |
| `scripts/` | Analysis scripts, numbered in run order |
| `results/` | Outputs written by scripts 23–29 |
| `tables/` | Source files for Tables 1–3, Supporting Information tables and Figure 1 (script 30) |
| `figures/` | Figure 1 in SVG, PDF, TIFF and PNG (script 31) |

## Scripts

| Script | Analysis |
|---|---|
| `23_rerun_frozen_no_proxy_uvmr.R` | Primary univariable MR (IVW, weighted median, MR-Egger, simple and weighted mode), Cochran's Q, MR-Egger intercept, MR-PRESSO, leave-one-out, single-SNP and Steiger analyses for the 10 outcomes, using directly matched SNPs only |
| `24_negative_control_analysis.R` | Negative-control outcomes (ease of skin tanning, `ukb-b-533`; natural blonde hair colour, `ukb-d-1747_1`) |
| `25_negative_control_influence_diagnostics.R` | Influence diagnostics for the negative-control analyses |
| `26_smoking_sensitivity_feasibility.R` | Instrument construction and feasibility screening for the multivariable MR (MVMR) models |
| `27_smoking_sensitivity_stage_b.R` | MVMR estimates and conditional F statistics (stair climbing + BMI + one smoking phenotype) |
| `28_smoking_covariance_sensitivity.R` | Conditional F statistics under source-informed sampling-error covariance scenarios |
| `29_targeted_pleiotropy_sensitivity.R` | Exclusion of adiposity- or education-associated instruments, exact same-size exclusion analysis and instrument influence summary |
| `30_build_bmc_main_table_sources.R` | Assembles table and figure source files from the results |
| `31_plot_figure1.R` | Draws Figure 1 |

## Running the analysis

Run the scripts in numerical order from the repository root:

```bash
for f in scripts/2[3-9]_*.R scripts/3[01]_*.R; do Rscript "$f"; done
```

All scripts run offline from `inputs/`. Scripts 24, 26 and 27 can re-download their inputs from the IEU OpenGWAS API when called with `--refresh`; this requires an OpenGWAS token saved in `~/.config/opengwas_token` and replaces the frozen inputs, so it is not needed to reproduce the reported results.

Random-number seeds for MR-PRESSO are fixed within the scripts. Scripts 23–30 write a run log and a manifest with SHA-256 checksums of their inputs and outputs to `results/` or `tables/`. MR-Egger and mode-based confidence intervals use t distributions matching their P values.

## Software

Analyses were run in R 4.5.2 with TwoSampleMR 0.6.22, MRPRESSO 1.0, MVMR 0.4.6, ieugwasr 1.1.0 and digest. Figure 1 additionally uses dplyr, ggplot2, patchwork, readr, scales, svglite and ragg.

## Data sources

| Trait | OpenGWAS ID | Source |
|---|---|---|
| Stair-climbing frequency (exposure) | `ukb-b-298` | UK Biobank |
| COPD differential diagnosis | `ukb-d-COPD_EXCL` | UK Biobank |
| COPD hospitalization; early-onset, late-onset and general COPD; emphysema; COPD-related respiratory insufficiency | `finn-b-COPD_HOSPITAL`, `finn-b-COPD_EARLY`, `finn-b-COPD_LATER`, `finn-b-J10_COPD`, `finn-b-J10_EMPHYSEMA`, `finn-b-COPD_INSUFFICIENCY` | FinnGen |
| FEV1, FVC, PEF | `ebi-a-GCST007432`, `ebi-a-GCST007429`, `ebi-a-GCST007430` | UK Biobank / GWAS Catalog |
| BMI (MVMR) | `ieu-a-835` | GIANT 2015 |
| Cigarettes per day; smoking initiation (MVMR) | `ieu-b-25`, `ieu-b-4877` | GSCAN |
| Cigarettes per day (MVMR) | `ieu-a-961` | TAG |

Use of these summary statistics is subject to the terms of the original data providers.

`inputs/instrument_annotation.csv` summarizes cataloged trait associations of the 16 instruments (NHGRI-EBI GWAS Catalog, IEU OpenGWAS and Ensembl, retrieved 8 August 2026). Its `判定` column holds the curator's classification label and `in_paper_16` marks the 16 analysed instruments (`是` = yes).

## License

Code is released under the MIT License (see `LICENSE`).
