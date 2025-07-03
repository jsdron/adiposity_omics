#########################################
# Script: dendo-heat_figure.R
# Description: Generates heatmaps for significant adiposity–analyte associations 
#              from both proteomic and metabolomic datasets.
# Key Outputs:
#   - Supplemental Figure 3 (heatmaps)
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(grid)
library(extrafont)

# Load fonts for PDF output
loadfonts(device = "pdf")

###### LOAD DATAFRAME(S) ######
# Load linear model results
prot <- fread(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.20250127.tsv.gz"))
prot$omic <- "Proteomic"

metab <- fread(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.20250127.tsv.gz"),
               select = c("term", "estimate", "std.error", "statistic", "p.value", "outcome", "adi_label"))
metab$omic <- "Metabolomic"

# Combine results
output_all <- rbind(prot, metab)

# Threshold and significance filter
threshold <- 41 + 445
output_all$sig <- ifelse(output_all$p.value < 0.05 / threshold, "sig", "ns")

# Standardize adiposity trait labels
output_all$adi_label <- ifelse(output_all$outcome == "vatadjbmi3", "VAT",
                         ifelse(output_all$outcome == "asatadjbmi3", "ASAT", "GFAT"))

###### METABOLOMIC HEATMAP DATA PREP ######
# Add metabolite annotations
metab_anno <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/input/metab_anno.txt", 
                    header = FALSE, stringsAsFactors = FALSE)
colnames(metab_anno) <- c("id", "outcome", "class", "unit", "lipo_size", "lipo_frac", 
                          "lipid_type", "general_type", "lipid_groups", "lipo_groups", 
                          "detail_groups", "label_name")
metab_anno$outcome <- str_replace(metab_anno$outcome, "met-d-IDL_IDL", "met-d-IDL")
metab_anno$outcome <- str_replace(metab_anno$outcome, "met-d-", "")
metab_anno <- metab_anno[metab_anno$unit != "ratio", ]
metab_anno <- metab_anno[, -1]  # Drop ID column

# Merge annotations
metab <- merge(metab, metab_anno, by.x = "term", by.y = "outcome")

# Create matrix: rows = analytes, columns = depots
matrix_metab <- metab %>%
  select(term, adi_label, estimate, general_type) %>%
  pivot_wider(names_from = adi_label, values_from = estimate) %>%
  arrange(general_type)

rownames(matrix_metab) <- matrix_metab$term
matrix_metab <- matrix_metab[, c("ASAT", "GFAT", "VAT")]  # Select columns in preferred order
matrix_metab <- as.matrix(matrix_metab)

# Generate metabolomic heatmap (rows clustered)
heat_metab <- pheatmap(matrix_metab,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       scale = "none",
                       legend = TRUE,
                       treeheight_row = 75,
                       treeheight_col = 0,
                       color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
                       fontsize = 12,
                       fontsize_col = 12,
                       fontsize_row = 12,
                       show_rownames = FALSE,
                       show_colnames = TRUE)

###### PROTEOMIC HEATMAP DATA PREP ######
matrix_prot <- prot %>%
  select(term, adi_label, estimate) %>%
  pivot_wider(names_from = adi_label, values_from = estimate)

rownames(matrix_prot) <- matrix_prot$term
matrix_prot <- matrix_prot[, c("ASAT", "GFAT", "VAT")]
matrix_prot <- as.matrix(matrix_prot)

# Generate proteomic heatmap (rows clustered)
heat_prot <- pheatmap(matrix_prot,
                      cluster_rows = TRUE,
                      cluster_cols = FALSE,
                      scale = "none",
                      legend = TRUE,
                      treeheight_row = 75,
                      treeheight_col = 0,
                      color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
                      fontsize = 12,
                      fontsize_col = 12,
                      fontsize_row = 12,
                      show_rownames = FALSE,
                      show_colnames = TRUE)

###### COMBINE PANELS AND EXPORT ######
heat_metab_grob <- arrangeGrob(heat_metab$gtable,
                               top = textGrob("A", x = unit(0, "npc"), just = "left",
                                              gp = gpar(fontsize = 14, fontface = "bold")))

heat_prot_grob <- arrangeGrob(heat_prot$gtable,
                              top = textGrob("B", x = unit(0, "npc"), just = "left",
                                             gp = gpar(fontsize = 14, fontface = "bold")))

# Save to PDF
pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure3.pdf"),
    width = 180 / 25.4, height = 5, family = "Arial")
grid.arrange(heat_metab_grob, heat_prot_grob, ncol = 2)
dev.off()
