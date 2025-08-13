#########################################
# Script: epi_linearReg_prot.R
# Description: Performs linear regression analyses (all sexes and sex-stratified) 
#              between proteins and MRI-derived adiposity traits 
#              (adjusted for BMI and height), and generates figures.
# Key Outputs: 
#   - Linear model results for protein associations
#   - Figure (protein Miami plot)
#   - Linear model results for sex-stratified protein associations
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
library(ggtext)
library(ggrepel)

# Load fonts for PDF output
loadfonts(device = "pdf")  

###### LOAD DATAFRAME(S) ######
# Load baseline dataset
baseline_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.noImpute.n5023.20240403.tsv.gz') 
colnames(baseline_prot) <- gsub("-", "_", colnames(baseline_prot))

###### SET VARIABLES ######
# Date the analysis was completed 
date <- "20250127"

# Fat depot variable names
adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

# Create the list of analyte variables
analyte <- colnames(baseline_prot)[which(colnames(baseline_prot) == "ERP44"):which(colnames(baseline_prot) == "COMMD1")]

# Correction factor
threshold <- 41+445 

###### FEATURE ASSOCIATION ######
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], 
                                    " + sex + age_instance2 + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_prot)
    tmp <- tidy(model)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}

# Label adiposity trait
output_prot$adi_label <- ifelse(output_prot$outcome == "vatadjbmi3", "VAT",
                                 ifelse(output_prot$outcome == "asatadjbmi3", "ASAT", "GFAT"))

###### SAVE RESULTS ######
write.table(output_prot, 
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.",date,".tsv.gz")), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

###### VISUALIZATION ######
# Load protein annotation file (chromosome and gene start position)
link <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/olink_protein_map_3k/olink_protein_map_3k_v1.tsv", 
              select = c("Assay", "chr", "gene_start"))
link$Assay <- gsub("-", "_", link$Assay)
link <- distinct(link, Assay, .keep_all = TRUE)

# Merge annotation info
output_prot <- merge(output_prot, link, by.x = "term", by.y = "Assay", all.x = TRUE)

# Format chromosome factor
output_prot$chr <- factor(output_prot$chr, levels = c(1:22, "X"))

# Signed -log10(p) value: negative for positive effects
output_prot$p_min <- ifelse(output_prot$estimate > 0, 
                            -log10(output_prot$p.value), 
                             log10(output_prot$p.value))

# Calculate cumulative positions for plotting
data_cum <- output_prot %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(gene_start, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  dplyr::select(chr, bp_add)

output_prot <- output_prot %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(bp_cum = gene_start + bp_add)

# Get chromosome label positions
axis_set <- output_prot %>%
  group_by(chr) %>%
  summarise(center = mean(bp_cum))

# Significance threshold
sig <- 0.05 / threshold

# Set y-axis limit for plotting
ylim <- output_prot %>%
  filter(p.value == min(p.value)) %>%
  mutate(ylim = abs(floor(log10(p.value))) + 2) %>%
  pull(ylim)

### Plot for ASAT
df_sub <- output_prot[output_prot$outcome == "asatadjbmi3", ]
df_sub$col <- factor(ifelse(df_sub$p.value < sig &
                              df_sub$term %in% output_prot[output_prot$outcome %in% c("gfatadjbmi3", "vatadjbmi3") & 
                                                           output_prot$p.value < sig, ]$term,
                            "A", ifelse(df_sub$p.value < sig, "B", df_sub$chr)),
                     levels = c(1:23, "A", "B"))
df_sub$Protein_lab <- ifelse(df_sub$p.value < 5e-15, df_sub$term, "")

plot_asat <- ggplot(df_sub, aes(x = bp_cum, y = p_min, color = col)) +
  geom_hline(yintercept = c(-log10(sig), log10(sig)), color = "grey40") + 
  geom_point(alpha = 1, size = 1.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = c(rep(c("grey90", "grey80"), 11), "grey90", "black", "#ECB41F")) +
  labs(x = "Genomic position", y = "-log<sub>10</sub>(P)") + 
  theme_cowplot() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = ggtext::element_markdown(size = 7),
    axis.title.x = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  ) +
  ggrepel::geom_text_repel(aes(label = Protein_lab), size = 5 / .pt)  # 5 pt converted to mm

### Plot for GFAT
df_sub <- output_prot[output_prot$outcome == "gfatadjbmi3", ]
df_sub$col <- factor(ifelse(df_sub$p.value < sig &
                              df_sub$term %in% output_prot[output_prot$outcome %in% c("asatadjbmi3", "vatadjbmi3") & 
                                                           output_prot$p.value < sig, ]$term,
                            "A", ifelse(df_sub$p.value < sig, "B", df_sub$chr)),
                     levels = c(1:23, "A", "B"))
df_sub$Protein_lab <- ifelse(df_sub$p.value < 5e-15, df_sub$term, "")

plot_gfat <- ggplot(df_sub, aes(x = bp_cum, y = p_min, color = col)) +
  geom_hline(yintercept = c(-log10(sig), log10(sig)), color = "grey40") + 
  geom_point(alpha = 1, size = 1.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = c(rep(c("grey90", "grey80"), 11), "grey90", "black", "#2883B1")) +
  labs(x = "Genomic position", y = "-log<sub>10</sub>(P)") + 
  theme_cowplot() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = ggtext::element_markdown(size = 7),
    axis.title.x = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  ) +
  ggrepel::geom_text_repel(aes(label = Protein_lab), size = 5 / .pt)  # 5 pt converted to mm

### Plot for VAT
df_sub <- output_prot[output_prot$outcome == "vatadjbmi3", ]
df_sub$col <- factor(ifelse(df_sub$p.value < sig &
                              df_sub$term %in% output_prot[output_prot$outcome %in% c("asatadjbmi3", "gfatadjbmi3") & 
                                                           output_prot$p.value < sig, ]$term,
                            "A", ifelse(df_sub$p.value < sig, "B", df_sub$chr)),
                     levels = c(1:23, "A", "B"))
df_sub$Protein_lab <- ifelse(df_sub$p.value < 5e-15, df_sub$term, "")

plot_vat <- ggplot(df_sub, aes(x = bp_cum, y = p_min, color = col)) +
  geom_hline(yintercept = c(-log10(sig), log10(sig)), color = "grey40") + 
  geom_point(alpha = 1, size = 1.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = c(rep(c("grey90", "grey80"), 11), "grey90", "black", "#C84D4C")) +
  labs(x = "Genomic position", y = "-log<sub>10</sub>(P)", title = "VATadjBMI") + 
  theme_cowplot() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = ggtext::element_markdown(size = 7),
    axis.title.x = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  ) +
  ggrepel::geom_text_repel(aes(label = Protein_lab), size = 5 / .pt)  # 5 pt converted to mm

### Save all panels in a single PDF
#pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure2.pdf"), 
#    width = 180 / 25.4, height = 220 / 25.4, family = "Arial")
plot_asat <- plot_asat + 
  ggtitle("ASAT") +
  theme(plot.title = element_text(size = 7))

plot_gfat <- plot_gfat + 
  ggtitle("GFAT") +
  theme(plot.title = element_text(size = 7))

plot_vat <- plot_vat + 
  ggtitle("VAT") +
  theme(plot.title = element_text(size = 7))

svg(
  paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure2.svg"),
  width  = 180 / 25.4,  
  height = 5,
  family = "Arial"
)

ggarrange(
  plot_asat,
  plot_gfat,
  plot_vat,
  labels = c("A", "B", "C"),
  ncol = 1, nrow = 3,
  font.label = list(size = 7, face = "bold", color = "black")
)

dev.off()

###### SEX-STRATIFIED FEATURE ASSOCIATION ######
output_prot_sex <- data.frame()

for (sex_group in c("F", "M")) {
  strat_data <- subset(baseline_prot, sex == sex_group)
  
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
      output_prot_sex <- rbind(output_prot_sex, result_row)
    }
  }
}

# Label adiposity trait
output_prot_sex$adi_label <- ifelse(output_prot_sex$outcome == "vatadjbmi3", "VAT",
                                     ifelse(output_prot_sex$outcome == "asatadjbmi3", "ASAT", "GFAT"))

###### SAVE SEX-STRATIFIED RESULTS ######
write.table(output_prot_sex,
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.sex_stratified.n5023.", date, ".tsv.gz")),
            sep = "\t", col.names = TRUE, row.names = FALSE)