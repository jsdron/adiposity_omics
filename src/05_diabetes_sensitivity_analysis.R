#########################################
# Script: diabetes_sensitivity_analysis.R
# Description: Sensitivity analysis excluding individuals with pre-diabetes or diabetes.
#              Compares analyte associations with adiposity traits against primary results.
# Key Outputs:
#   - Supplemental Figure 6 (scatter plots of effect estimates)
#   - Supplemental Table (analyte comparison results)
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

font_import()  # This imports all available fonts to R
loadfonts(device = "pdf")  # Load fonts for PDF outputt

###### PREPARE DATAFRAME(S) ######
# Load baseline dataset
baseline_prot <- fread("/Users/michael.tian/Desktop/Natarajan_Lab/adiposity/adiposity_omics/data/baseline_adi_prot.noImpute.n5023.20240730.tsv.gz")
colnames(baseline_prot) <- gsub("-", "_", colnames(baseline_prot))

baseline_metab <- fread("/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/data/baseline_adi_metab.n22630.20240730.tsv.gz")

# Adding labels for pre-diabetes and diabetes based on HbA1C
hba1c <- fread("/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/data/UKBB_TG_BackgroundVariables_22020.txt", select = c('eid', 'HgbA1c'))

baseline_prot <- merge(baseline_prot, hba1c, by = "eid")
baseline_metab <- merge(baseline_metab, hba1c, by = "eid")

# Add diabetes labels
baseline_prot$Diabetes <- ifelse(baseline_prot$HgbA1c < 38.8, "Normal", 
                          ifelse(baseline_prot$HgbA1c < 47.5, "Pre-Diabetes", "Diabetes"))

baseline_metab$Diabetes <- ifelse(baseline_metab$HgbA1c < 38.8, "Normal", 
                                 ifelse(baseline_metab$HgbA1c < 47.5, "Pre-Diabetes", "Diabetes"))

# Filtering down to individuals without pre-diabetes or diabetes
baseline_noDM <- baseline_prot[baseline_prot$Diabetes == "Normal", ] 
baseline_metab <- baseline_metab[baseline_metab$Diabetes == "Normal", ] 

###### SET VARIABLES ######
# Date the analysis was completed 
date <- "20250127"

# Fat depot variable names
adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

# Create the list of analyte variables
analyte_prot <- colnames(baseline_prot)[which(colnames(baseline_prot) == "ERP44"):which(colnames(baseline_prot) == "COMMD1")]

analyte <- colnames(baseline_metab)[which(colnames(baseline_metab) == "Clinical_LDL_C"):which(colnames(baseline_metab) == "Omega_6_pct_PUFA")]

adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height


# Correction factor
threshold <- 41+445 

###### FEATURE ASSOCIATION ######
# Proteins
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte_prot)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte_prot[j], "+ sex + age_instance2 + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_noDM) 
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte_prot[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}

output_prot$adi_label <- ifelse(output_prot$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

output_prot$DMstatus <- "No diabetes"

# Metabolites
output_metab <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + age_instance2 + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_metab)
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_metab <- rbind(output_metab, result_row)
  }
}

output_metab$adi_label <- ifelse(output_metab$outcome == "vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_metab$outcome == "asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

# For grouping purposes, we need to add the "super groups" for each metabolite. _pct values removed
raw <- fread(file = "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/data/raw_jwl_jd.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(raw) = c("id", "outcome", "class", "unit", "lipo_size", "lipo_frac", "lipid_type", "general_type", "lipid_groups", "lipo_groups", "detail_groups", "label_name")
raw$outcome <- str_replace(raw$outcome, "met-d-IDL_IDL", "met-d-IDL")
raw$outcome <- str_replace(raw$outcome, "met-d-", "")
raw <- raw[, -c(1)]

# Exclude ratios and pct
raw <- raw[raw$unit!="ratio",]

met_list <- raw$outcome

output_metab <- merge(output_metab, raw, by.x = "term", by.y = "outcome")

rm(raw)

###### SAVE RESULTS ######
write.table(output_prot, 
            file = gzfile(paste0("/Users/michael.tian/Desktop/Natarajan_Lab/adiposity/adiposity_omics/results/lmResults_adi_prot_nodiabetes.n5023.",date,".tsv.gz")), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

write.table(output_metab,
            file = gzfile(paste0("/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/results/lmResults_adi_metab_nodiabetes.n22630.",date,".tsv.gz")), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


###### LOAD DATAFRAMES ######
prot_dm <- output_prot
prot_dm$label <- paste0(prot_dm$term, "_", prot_dm$outcome)
prot_main <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.20250127.tsv.gz')
prot_main$label <- paste0(prot_main$term, "_", prot_main$outcome)

epi_dm <- output_metab
epi_dm$label <- paste0(epi_dm$term, "_", epi_dm$outcome)
epi_main <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.20250127.tsv.gz')
epi_main$label <- paste0(epi_main$term, "_", epi_main$outcome)


###### COMPARISONS ######
metab <- merge(epi_main, epi_dm, by = "label")

metab$adi_label <- ifelse(metab$outcome.x == "vatadjbmi3", "VAT",
                          ifelse(metab$outcome.x == "asatadjbmi3", "ASAT", "GFAT"))

metab$compare <- ifelse(metab$p.value.x < 0.05/threshold & metab$p.value.y >= 0.05/threshold, "Primary only",
                        ifelse(metab$p.value.x >= 0.05/threshold & metab$p.value.y < 0.05/threshold, "Sensitivity only",
                              ifelse(metab$p.value.x < 0.05/threshold & metab$p.value.y < 0.05/threshold, "Both", "Not significant")))
metab$compare <- factor(metab$compare, levels = c("Not significant", "Both", "Primary only", "Sensitivity only"))

metab$sens_sig <- ifelse(metab$p.value.y < 0.05/threshold, "sig", "ns")

table(metab$adi_label, metab$sens_sig)

###### SUPPLEMENTAL FIGURE 6 ######
# Metabolite plot
metab_plot <- ggplot(metab, aes(x = estimate.x, y = estimate.y)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "#dddddd") +
  geom_vline(xintercept = 0, linetype = "solid", color = "#dddddd") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5, colour = "#dddddd") +
  geom_point(data = subset(metab, !(compare %in% c("Sensitivity only", "Primary only"))),
             aes(color = compare), size = 1, alpha = 0.7, show.legend = TRUE) +
  geom_smooth(method = "lm", formula = y ~ x, aes(group = 1), colour = "black", linewidth = 0.5, fullrange = TRUE) +
  geom_point(data = subset(metab, compare %in% c("Sensitivity only", "Primary only")),
             aes(color = compare), size = 1, alpha = 0.7, show.legend = TRUE) +
  scale_color_manual(values = c("Not significant" = "lightgray",
                                "Both" = "#693182",
                                "Sensitivity only" = "darkblue",
                                "Primary only" = "darkred")) +
  facet_wrap(~adi_label) +
  labs(x = "Primary effect estimate", y = "Sensitivity effect estimate", color = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.justification = "center") +
  ggpubr::stat_cor(r.accuracy = 0.01, cor.coef.name = "r", position = "identity",
                   aes(x = estimate.x, y = estimate.y,
                       label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                   inherit.aes = FALSE)

# Protein plot
prot <- merge(prot_main, prot_dm, by = "label")

prot$adi_label <- ifelse(prot$outcome.x == "vatadjbmi3", "VAT",
                  ifelse(prot$outcome.x == "asatadjbmi3", "ASAT", "GFAT"))

prot$compare <- ifelse(prot$p.value.x < 0.05/threshold & prot$p.value.y >= 0.05/threshold, "Primary only",
                ifelse(prot$p.value.x >= 0.05/threshold & prot$p.value.y < 0.05/threshold, "Sensitivity only",
                ifelse(prot$p.value.x < 0.05/threshold & prot$p.value.y < 0.05/threshold, "Both", "Not significant")))
prot$compare <- factor(prot$compare, levels = c("Not significant", "Both", "Primary only", "Sensitivity only"))

prot$sens_sig <- ifelse(prot$p.value.y < 0.05/threshold, "sig", "ns")

table(prot$adi_label, prot$sens_sig)

prot_plot <- ggplot(prot, aes(x = estimate.x, y = estimate.y)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "#dddddd") +
  geom_vline(xintercept = 0, linetype = "solid", color = "#dddddd") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5, colour = "#dddddd") +
  geom_point(data = subset(prot, !(compare %in% c("Sensitivity only", "Primary only"))),
             aes(color = compare), size = 1, alpha = 0.7, show.legend = TRUE) +
  geom_smooth(method = "lm", formula = y ~ x, aes(group = 1), colour = "black", linewidth = 0.5, fullrange = TRUE) +
  geom_point(data = subset(prot, compare %in% c("Sensitivity only", "Primary only")),
             aes(color = compare), size = 1, alpha = 0.7, show.legend = TRUE) +
  scale_color_manual(values = c("Not significant" = "lightgray",
                                "Both" = "#693182",
                                "Sensitivity only" = "darkblue",
                                "Primary only" = "darkred")) +
  facet_wrap(~adi_label) +
  labs(x = "Primary effect estimate", y = "Sensitivity effect estimate", color = "Significance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.justification = "center") +
  ggpubr::stat_cor(r.accuracy = 0.01, cor.coef.name = "r", position = "identity",
                   aes(x = estimate.x, y = estimate.y,
                       label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                   inherit.aes = FALSE)

###### PLOT EXPORT ######
pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure4.pdf", 
    width = 180 / 25.4, height = 180 / 25.4, family = "Arial")

ggarrange(metab_plot, prot_plot,
          labels = c("A", "B"),
          legend = "bottom", 
          common.legend = TRUE,
          nrow = 2)

dev.off()

###### SUPPLEMENTAL TABLE ######
metab_cast <- dcast(metab, term.x ~ outcome.x, value.var = "compare")
metab_cast$group <- "Metabolomic"

prot_cast <- dcast(prot, term.x ~ outcome.x, value.var = "compare")
prot_cast$group <- "Proteomic"

merged_cast <- rbind(metab_cast, prot_cast)
colnames(merged_cast) <- c("Analyte", "ASAT", "GFAT", "VAT", "Omic")
merged_cast <- merged_cast[, c("Omic", "Analyte", "ASAT", "GFAT", "VAT")]


###### SAVE RESULTS ######
write.table(merged_cast,
            file = gzfile("/Users/michael.tian/Desktop/Natarajan_Lab/adiposity/adiposity_omics/results/supplemental_table_diabetes_sensitivity.tsv.gz"),
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
