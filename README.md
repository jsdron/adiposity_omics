# Integrative Multi-Omic Analysis of Regional Fat Depots

This repository contains analysis scripts accompanying the manuscript:

**Dron et al.** _Integrative multi-omic analysis of regional fat depots reveals molecular determinants of human fat distribution_ (2025). [Submitted].

## Overview

We integrated abdominal MRI-derived fat depot volumes (VAT, ASAT, GFAT) with circulating metabolomics and proteomics from UK Biobank to identify molecular signatures of regional adiposity. Using association testing, sex-stratified and diabetes-free sensitivity analyses, pathway enrichment, and bidirectional Mendelian randomization (MR), we nominate candidate circulating factors that influence fat distribution and are linked to cardiometabolic disease risk.

Note that access to the UK Biobank is required to replicate these results.

## Directory Structure

```bash
./src/                      # Analysis scripts (see detailed list below)
./sumstats/                 # GWAS summary statistics
```

## Script Index
| Script                                  | Description                                                                                                                                           |
| --------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01_prepare_baseline_datasets.R`        | Prepares baseline data by filtering samples, adjusting adiposity traits, imputing and scaling omics variables, and saving cleaned datasets.           |
| `02_variance_explained.R`               | Computes the number of principal components needed to explain 90% and 99% of variance in proteomics data.                                             |
| `03a_epi_linearReg_metab.R`             | Performs linear regression of metabolites with VAT, ASAT, and GFAT (adjusted for BMI and height), including sex-stratified models and visualizations. |
| `03b_epi_linearReg-interaction_metab.R` | Tests for sex interactions between metabolites and fat depots; compares formal interaction models with sex-stratified results.                        |
| `04a_epi_linearReg_prot.R`              | Performs linear regression of proteins with VAT, ASAT, and GFAT (adjusted for BMI and height), including sex-stratified models and visualizations.    |
| `04b_epi_linearReg-interaction_prot.R`  | Tests for sex interactions between proteins and fat depots; compares formal interaction models with sex-stratified results.                           |
| `05_diabetes_sensitivity_analysis.R`    | Repeats metabolite and protein association analyses after excluding individuals with pre-diabetes or diabetes; compares with primary results.         |
| `06_table1.R`                           | Generates summary Table 1 of clinical characteristics by omics group (metabolomics vs. proteomics), including diabetes and WHR status.                |
| `07_summary_figures.R`                  | Summarizes significant analyte-depot associations using bar plots and upset plots to visualize shared and unique signals.                             |
| `08_dendo-heat_figure.R`                | Generates heatmaps of standardized effect estimates for significant analyte-depot associations.                                                       |
| `09_sex-strat_figure.R`                 | Generates volcano and scatter plots comparing male and female analyte associations across depots.                                                     |
| `10_pathway_figure.R`                   | Generates pathway enrichment plots for proteins associated with ASAT, VAT, GFAT, and their overlaps using g\:Profiler.                                |
| `11_Out-Fat_Exp-Metab.R`                | Runs two-sample MR of metabolite exposures on fat depots.                                                                                             |
| `12_Out-Fat_Exp-Prot.R`                 | Runs two-sample MR of protein exposures on fat depots.                                                                                                |
| `13_Out-Metab_Exp-Fat.R`                | Runs two-sample MR of fat depots on metabolite levels.                                                                                                |
| `14_Out-Prot_Exp-Fat.R`                 | Runs two-sample MR of fat depots on protein levels.                                                                                                   |
| `15_MR_outputs.R`                       | Organizes and formats MR results across runs; harmonizes instruments and effect directions.                                                           |
| `16_MR_figure.R`                        | Generates forest plots for forward and reverse MR results by fat depot.                                                                               |
| `17_coxph_prot_candidates.R`            | Fits Cox models for protein MR candidates and plots associations with coronary artery disease and type 2 diabetes.                                    |
| `18_coxph_metab_candidates.R`           | Fits Cox models for metabolite MR candidates and plots associations with coronary artery disease and type 2 diabetes.                                 |

## Citation

If you use this code, please cite:

> Dron JS, Pan M, Schuermans A, et al. Integrative multi-omic analysis of regional fat depots reveals molecular determinants of human fat distribution. 2025 (Submitted).
