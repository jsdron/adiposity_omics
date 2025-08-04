#########################################
# Script: epi_linearReg-interaction_metab.R
# Description: Tests for sex interactions between MRI-derived adiposity 
#              traits and metabolites. Compares formal interaction 
#              results with stratified models, and summarizes patterns by sex.
# Key Outputs:
#   - Linear model results with interaction terms
#   - Interaction summary merged with stratified models
#   - Summary counts and percentages by sex and fat depot
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tableone)
library(stringr)
library(broom)
library(ggpubr)
library(extrafont)

# Imports all available fonts to R
font_import()  

# Load fonts for PDF output
loadfonts(device = "pdf")  

###### LOAD DATAFRAME(S) ######
# Load baseline dataset
baseline_metab <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.n22630.20240403.tsv.gz")

###### SET VARIABLES ######
# Date the analysis was completed
date <- "20250411"

# Fat depot variable names
adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3")

# Create the list of analyte variables
analyte <- colnames(baseline_metab)[which(colnames(baseline_metab) == "Clinical_LDL_C"):which(colnames(baseline_metab) == "Omega_6_pct_PUFA")]

# Correction factor for multiple testing
threshold <- 41 + 445

###### FEATURE ASSOCIATION (INTERACTION MODELS) ######
output_metab <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(
      adi_traits[i], " ~ ", analyte[j],
      " + sex + age_instance2 + mriNum + time_between + sex*", analyte[j])
    )
    model <- lm(lm_formula, data = baseline_metab)
    tmp <- tidy(model)
    result_row <- tmp[which(tmp$term == paste0(analyte[j], ":sexM")), ]
    result_row$outcome <- adi_traits[i]
    result_row$analyte <- analyte[j]
    output_metab <- rbind(output_metab, result_row)
  }
}

# Label adiposity trait
output_metab$adi_label <- ifelse(output_metab$outcome == "vatadjbmi3", "VAT",
                                 ifelse(output_metab$outcome == "asatadjbmi3", "ASAT", "GFAT"))

# Identify significant interactions
output_metab$sig_interaction <- ifelse(output_metab$p.value < 0.05 / threshold, "Significant", "NS")

###### ADD ANNOTATIONS ######
raw <- fread("/Volumes/medpop_esp2/jdron/projects/cihr_metab/analysis/v1_preNov2023/02_MR/data/raw_jwl_jd.txt", 
             header = FALSE, stringsAsFactors = FALSE)
colnames(raw) <- c(
  "id", "outcome", "class", "unit", "lipo_size", "lipo_frac", 
  "lipid_type", "general_type", "lipid_groups", "lipo_groups", 
  "detail_groups", "label_name")
raw$outcome <- str_replace(raw$outcome, "met-d-IDL_IDL", "met-d-IDL")
raw$outcome <- str_replace(raw$outcome, "met-d-", "")
raw <- raw[, -1]

# We are not considering the ratio measures
raw <- raw[raw$unit != "ratio", ]

# Merge annotation info
output_metab <- merge(output_metab, raw, by.x = "term", by.y = "outcome")

# Identify significant interactions
output_metab$sig_interaction <- ifelse(output_metab$p.value <0.05/threshold, "Significant", "NS")

###### SAVE RESULTS ######
write.table(output_metab,
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.interaction.", date, ".tsv.gz")),
            sep = "\t", col.names = TRUE, row.names = FALSE)

output_metab$Stratification.Group <- "All"

###### INTERSECT WITH SEX-STRATIFIED RESULTS ######
# Load interaction results
int_met <- output_metab
int_met$label <- paste0(int_met$adi_label, "_", int_met$analyte)
int_met <- int_met[, c("label", "p.value", "sig_interaction")]
colnames(int_met)[2] <- "Interaction P-value"

# Load stratified results
sex_strat_met <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.sex_stratified.n22630.20240403.tsv.gz")
sex_strat_met$label <- paste0(sex_strat_met$adi_label, "_", sex_strat_met$analyte)

# Merge interaction and stratified data
test_met <- merge(sex_strat_met, int_met, by = "label", all = TRUE)

# Identify metabolites wiht significant sex* interaction
sig_interaction_met <- subset(test_met, test_met$sig_interaction == "Significant")
sig_interaction_met <- subset(sig_interaction_met, sig_interaction_met$Stratification.Group != "All")

# Get number of unique sex-interacting metabolites per depot
total_by_depot <- sig_interaction_met %>%
  filter(Stratification.Group %in% c("Male", "Female")) %>%
  group_by(outcome) %>%
  summarise(n_labels = n_distinct(label))

# Determine stronger sex effect
stronger_effects <- sig_interaction_met %>%
  filter(Stratification.Group %in% c("Male", "Female")) %>%
  mutate(abs_beta = abs(estimate)) %>%
  group_by(label, outcome) %>%
  filter(n() == 2) %>%
  summarise(
    stronger_in = case_when(
      abs_beta[Stratification.Group == "Male"] > abs_beta[Stratification.Group == "Female"] ~ "Male",
      abs_beta[Stratification.Group == "Male"] < abs_beta[Stratification.Group == "Female"] ~ "Female",
      TRUE ~ "Equal"
    ),
    .groups = "drop"
  )

# Count by depot
counts_by_depot <- stronger_effects %>%
  count(outcome, stronger_in) %>%
  pivot_wider(names_from = stronger_in, values_from = n, values_fill = 0)

# Calculate proportions
percent_by_depot <- counts_by_depot %>%
  left_join(total_by_depot, by = "outcome") %>%
  mutate(
    Male_percent = round(100 * Male / n_labels, 1),
    Female_percent = round(100 * Female / n_labels, 1)
  )

