#########################################
# Script: epi_linearReg_metab.R
# Description: Performs linear regression analyses (for all sexes and sex-stratified) between  
#              metabolites and MRI-derived adiposity traits 
#              (adjusted for BMI and height), adds metabolite annotations,
#              and generates figures.
# Key Outputs: 
#   - Linear model results for metabolite associaitons
#   - Figures (Supplemental Figure 4)
#   - Linear model results for sex-stratified metabolite associations
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tableone)
library(stringr)
library(cowplot)
library(broom)
library(ggpubr)
library(extrafont)

# Imports all available fonts to R
font_import(prompt=FALSE)  

# Load fonts for PDF output
loadfonts(device = "pdf")  

###### FUNCTIONS ######
# The function to generate forest plots for the metabolite associations
plot_forest <- function(df, ncol_wrap, title_input) {
  ggplot(data = df,
         aes(x = label_name, y = estimate, color = adi_label,
             fill = ifelse(p.value < 0.05/threshold, adi_label, "white"),
             shape = ifelse(p.value < 0.05/threshold, "Significant", "Not significant"))) +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(stroke = 1, show.legend = TRUE,
               position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "", y = "Effect estimate (SE)", color = "", title = title_input, shape = "", fill = "") +
    theme_bw() +
    coord_flip() +
    scale_shape_manual(values = c("Significant" = 16, "Not significant" = 21)) +
    theme(legend.position = "top", 
          title = element_text(size = 14, face = "bold", hjust=0),
          strip.text = element_text(size = 12),
          axis.title.x = element_text(size = 12, face="plain", hjust = 0.5),
          axis.text = element_text(size = 11), 
          legend.justification = "center") +
    scale_color_manual(values = c("#ECB41F", "#2883B1", "#C84D4C")) +
    scale_fill_manual(values = c("#ECB41F", "#2883B1", "#C84D4C", "white"), guide = "none") +
    facet_wrap(~detail_groups, scales = "free", ncol = ncol_wrap)
}

###### LOAD DATAFRAME(S) ######
# Load baseline dataset
baseline_metab <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.n22630.20240403.tsv.gz")

###### SET VARIABLES ######
# Date the analysis was completed 
date <- "20250127"

# Fat depot variable names
adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

# Create the list of analyte variables
analyte <- colnames(baseline_metab)[which(colnames(baseline_metab) == "Clinical_LDL_C"):which(colnames(baseline_metab) == "Omega_6_pct_PUFA")]

# Correction factor
threshold <- 41 + 445

###### FEATURE ASSOCIATION ######
output_metab <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], 
                                    "+ sex + age_instance2 + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_metab)
    tmp <- tidy(model)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_metab <- rbind(output_metab, result_row)
  }
}

# Label adiposity trait
output_metab$adi_label <- ifelse(output_metab$outcome == "vatadjbmi3", "VAT",
                                 ifelse(output_metab$outcome == "asatadjbmi3", "ASAT", "GFAT"))

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

###### SAVE RESULTS ######
write.table(output_metab,
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.", date, ".tsv.gz")),
            sep = "\t", col.names = TRUE, row.names = FALSE)

###### VISUALIZATION ######
output_metab <- subset(output_metab, term != "Unsaturation")

# Recode label names for clarity in the plot
output_metab$label_name <- recode(output_metab$label_name,
                                  "Sphingomyelins" = "SM",
                                  "Phosphoglyc" = "PhGlyc",
                                  "Phosphatidylc" = "PC",
                                  "Cholines" = "Cho",
                                  "Total_BCAA" = "BCAA",
                                  "Acetoacetate" = "AcAc",
                                  "bOHbutyrate" = "BHB",
                                  "Omega_3" = "Omega 3",
                                  "Omega_6" = "Omega 6"
                                  )

# Set factor levels
output_metab <- output_metab %>%
  arrange(detail_groups, desc(outcome)) %>%
  mutate(label_name = factor(label_name, levels = unique(label_name)),
         detail_groups = factor(detail_groups))

# GCreate the metabolite groups for plotting
a_groups <- c("Cholesterol","Cholesteryl esters","Free cholesterol",
              "Triglycerides","Phospholipids","Total lipids",
              "Particle concentration","Particle size","Apolipoproteins")
b_groups <- c("HDL - XLarge", "HDL - Large","HDL - Medium","HDL - Small",
              "IDL", "LDL - Large", "LDL - Medium","LDL - Small")
c_groups <- c("Chylomicron","VLDL - XLarge","VLDL - Large",
              "VLDL - Medium","VLDL - Small","VLDL - XSmall")
d_groups <- c('Amino acids','Fatty acids','Glycolysis')
e_groups <- c('Inflammation','Fluid balance','Ketone bodies')

# Generate the different plots based on different metabolite groups
a_panel <- plot_forest(subset(output_metab, output_metab$detail_groups %in% a_groups), 3, "Lipids")
b_panel <- plot_forest(subset(output_metab, output_metab$detail_groups %in% b_groups), 4, "Lipoproteins")
c_panel <- plot_forest(subset(output_metab, output_metab$detail_groups %in% c_groups), 3, NA)
d_panel <- plot_forest(subset(output_metab, output_metab$detail_groups %in% d_groups), 3, "Other metabolites")
e_panel <- plot_forest(subset(output_metab, output_metab$detail_groups %in% e_groups), 3, NA)

# Creating the plots in separate files, to be combined later
pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure4a.pdf"), 
    width = 180/25.4, height = 180/25.4,family = "Arial")
print(a_panel + theme(legend.position = "none"))
dev.off()

pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure4b.pdf"), 
    width = 180/25.4, height = 220/25.4,family = "Arial") 
ggarrange(
  b_panel + theme(axis.title.x = element_blank(),
                  legend.position = "none"),
  c_panel + theme(plot.title = element_blank(),
                  legend.position = "none"),
  ncol = 1,
  legend = "none",
  heights = c(1, 1)
)
dev.off()

pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure4c.pdf"), 
    width = 180/25.4, height = 140/25.4,family = "Arial") # 
ggarrange(
  d_panel + theme(axis.title.x = element_blank()),
  e_panel + theme(plot.title = element_blank()),
  ncol = 1,
  common.legend = TRUE,
  legend = "bottom",
  heights = c(2, 1)
)
dev.off()

###### SEX-STRATIFIED FEATURE ASSOCIATION ######
output_metab_sex <- data.frame()

for (sex_group in c("F", "M")) {
  strat_data <- subset(baseline_metab, sex == sex_group)
  
  for (i in 1:length(adi_traits)) {
    for (j in 1:length(analyte)) {
      lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j],
                                      " + age_instance2 + mriNum + time_between"))
      model <- lm(lm_formula, data = strat_data)
      tmp <- tidy(model)
      result_row <- tmp[which(tmp$term == analyte[j]), ]
      result_row$outcome <- adi_traits[i]
      result_row$analyte <- analyte[j]
      result_row$Stratification.Group <- ifelse(sex_group == "F", "Female", "Male")
      output_metab_sex <- rbind(output_metab_sex, result_row)
    }
  }
}
# Label adiposity trait
output_metab_sex$adi_label <- ifelse(output_metab_sex$outcome == "vatadjbmi3", "VAT",
                                 ifelse(output_metab_sex$outcome == "asatadjbmi3", "ASAT", "GFAT"))

###### ADD ANNOTATIONS ######
# Merge annotation info
output_metab_sex <- merge(output_metab_sex, raw, by.x = "term", by.y = "outcome")

###### SAVE RESULTS ######
write.table(output_metab_sex,
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.sex_stratified.n22630.", date, ".tsv.gz")),
            sep = "\t", col.names = TRUE, row.names = FALSE)