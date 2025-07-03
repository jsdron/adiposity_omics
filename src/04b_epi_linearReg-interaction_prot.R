#########################################
# Script: epi_linearReg-interaction_prot.R
# Description: Tests for sex interactions between MRI-derived adiposity 
#              traits and proteins Compares formal interaction 
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
baseline_prot <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.noImpute.n5023.20240403.tsv.gz")
colnames(baseline_prot) <- gsub("-", "_", colnames(baseline_prot))

###### SET VARIABLES ######
# Date the analysis was completed 
date <- "20250127"

# Fat depot variable names
adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

# Create the list of analyte variables
analyte <- colnames(baseline_prot)[which(colnames(baseline_prot) == "ERP44"):which(colnames(baseline_prot) == "COMMD1")]

# Correction factor
threshold <- 41 + 445 

###### FEATURE ASSOCIATION ######
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j],
                                    " + sex + age_instance2 + mriNum + time_between + sex*", analyte[j]))
    model <- lm(lm_formula, data = baseline_prot)
    tmp <- tidy(model)
    result_row <- tmp[which(tmp$term == paste0(analyte[j], ":sexM")), ]
    result_row$outcome <- adi_traits[i]
    result_row$analyte <- analyte[j]
    output_prot <- rbind(output_prot, result_row)
  }
}

# Label adiposity trait
output_prot$adi_label <- ifelse(output_prot$outcome == "vatadjbmi3", "VAT",
                                 ifelse(output_prot$outcome == "asatadjbmi3", "ASAT", "GFAT"))

# Mark significant interactions
output_prot$sig_interaction <- ifelse(output_prot$p.value < 0.05 / threshold, "Significant", "NS")

###### SAVE RESULTS ######
write.table(output_prot, 
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.interaction.",date,".tsv.gz")), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

output_prot$Stratification.Group <- "All"

###### INTERSECT WITH SEX-STRATIFIED RESULTS ######
# Load interaction results
int_prot <- output_prot
int_prot$label <- paste0(int_prot$adi_label, "_", int_prot$analyte)
int_prot <- int_prot[, c("label", "p.value", "sig_interaction")]
colnames(int_prot)[2] <- c("Interaction P-value")

# Load stratified results
sex_strat_prot <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.sex_stratified.n5023.20250127.tsv.gz")
sex_strat_prot$label <- paste0(sex_strat_prot$adi_label, "_", sex_strat_prot$analyte)

# Merge interaction and stratified data
test <- merge(sex_strat_prot, int_prot, by="label", all=TRUE)

# Identify proteims wiht significant sex* interaction
sig_interaction_prot <- subset(test, test$sig_interaction=="Significant")
sig_interaction_prot <- subset(sig_interaction_prot, sig_interaction_prot$`Stratification group`!="All")
colnames(sig_interaction_prot)[2] <- "Stratification.group"

# Get number of unique sex-interacting metabolites per depot
total_by_depot_prot <- sig_interaction_prot %>%
  filter(Stratification.group %in% c("Male", "Female")) %>%
  group_by(outcome) %>%
  summarise(n_labels = n_distinct(label), .groups = "drop")

# Determine stronger effect by sex
stronger_effects_prot <- sig_interaction_prot %>%
  filter(Stratification.group %in% c("Male", "Female")) %>%
  mutate(abs_beta = abs(estimate)) %>%
  group_by(label, outcome) %>%
  filter(n() == 2) %>%  # Only keep complete Male/Female pairs
  summarise(
    stronger_in = case_when(
      abs_beta[Stratification.group == "Male"] > abs_beta[Stratification.group == "Female"] ~ "Male",
      abs_beta[Stratification.group == "Male"] < abs_beta[Stratification.group == "Female"] ~ "Female",
      TRUE ~ "Equal"
    ),
    .groups = "drop"
  )

# Count stronger sex effects per depot
counts_by_depot_prot <- stronger_effects_prot %>%
  count(outcome, stronger_in) %>%
  pivot_wider(names_from = stronger_in, values_from = n, values_fill = 0)

# Merge and calculate percentages
percent_by_depot_prot <- counts_by_depot_prot %>%
  left_join(total_by_depot_prot, by = "outcome") %>%
  mutate(
    Male_percent = round(100 * Male / n_labels, 1),
    Female_percent = round(100 * Female / n_labels, 1)
  )