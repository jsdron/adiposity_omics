#########################################
# Script: summary_figures.R
# Description: Combines metabolomic and proteomic results to:
#              (1) compare significant analytes across fat depots, 
#              (2) visualize shared vs. distinct signals via barplots and upset plots.
# Key Outputs:
#   - Figure 2 (upset plot)
#   - Supplemental Figure 2 (barplots)
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
library(tidyverse)
library(ggtext)
library(ggpubr)
library(gtools)
library(ggrepel)
library(ggbreak)
library(UpSetR) 
library(gridExtra)
library(extrafont)

# Load fonts for PDF output
font_import()  # This imports all available fonts to R
loadfonts(device = "pdf")

###### SET VARIABLES ######
# Correction factor
threshold <- 41+445 

###### PREPARE DATAFRAMES ######
a <- fread(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.20250127.tsv.gz"))
a$omic <- "Proteomic"

b <- fread(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.20250127.tsv.gz"),
           select = c("term", "estimate", "std.error", "statistic", "p.value", "outcome", "adi_label"))
b$omic <- "Metabolomic"

# Combine
output_all <- rbind(a, b)

# Annotate
output_all$sig <- ifelse(output_all$p.value < 0.05 / threshold, "sig", "ns")
output_all$adi_label <- ifelse(output_all$outcome == "vatadjbmi3", "VAT",
                         ifelse(output_all$outcome == "asatadjbmi3", "ASAT", "GFAT"))

rm(a, b)

# Subset significant results
metab <- subset(output_all, omic == "Metabolomic" & sig == "sig")
prot <- subset(output_all, omic == "Proteomic" & sig == "sig")

###### SUPPLEMENTAL FIGURE 2: BARPLOTS OF SIGNIFICANT ANALYTES ######
plot_data <- output_all %>%
  count(adi_label, omic, sig) %>%
  group_by(adi_label, omic) %>%
  mutate(percentage = n / sum(n) * 100,
         rounded_per = round(percentage, 1))

# Absolute count plot
fig1a_i <- ggplot(plot_data[plot_data$sig == "sig", ], aes(x = adi_label, y = n, fill = omic)) +
  geom_bar(stat = "identity") +
  labs(x = "Fat depot", y = "Significant analytes (n)", fill = "") +
  scale_fill_manual(values = c("#E1D2B2", "#A68A6D")) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11)) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))

# Percent bar plot
fig1a_ii <- ggplot(plot_data[plot_data$sig == "sig", ], aes(x = adi_label, y = rounded_per, fill = omic)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Fat depot", y = "Significant analytes (%)", fill = "") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  scale_fill_manual(values = c("#E1D2B2", "#A68A6D")) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11)) +
  geom_text(aes(label = scales::percent(rounded_per / 100)), 
            position = position_dodge(width = 0.9), vjust = -0.5)

# Combine for Supplemental Figure 2
sFig1cd <- ggarrange(fig1a_i, fig1a_ii, ncol = 2, common.legend = TRUE,
                     legend = "bottom", align = "h", labels = c("A", "B"))

# Save PDF
pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure2.pdf"), 
    width = 180 / 25.4, height = 3, family = "Arial")
sFig1cd
dev.off()

###### FIGURE 2: UPSET PLOT ######
sig_analyte <- output_all[output_all$sig == "sig", ]

# Get significant analytes by depot
A <- sig_analyte[sig_analyte$adi_label == "VAT", ]$term
B <- sig_analyte[sig_analyte$adi_label == "ASAT", ]$term
C <- sig_analyte[sig_analyte$adi_label == "GFAT", ]$term

upset_df <- data.frame(Analyte = sig_analyte$term, Omic = sig_analyte$omic)
upset_df$VAT <- ifelse(upset_df$Analyte %in% A, 1, 0)
upset_df$ASAT <- ifelse(upset_df$Analyte %in% B, 1, 0)
upset_df$GFAT <- ifelse(upset_df$Analyte %in% C, 1, 0)

# Create upset plot
fig1b <- upset(upset_df,
               sets = c("VAT", "ASAT", "GFAT"),
               order.by = "freq",
               mainbar.y.label = "Intersection size",
               sets.x.label = "Significant analytes",
               empty.intersections = "on",
               query.legend = "bottom",
               text.scale = rep(1.5, 6),
               point.size = 3,
               sets.bar.color = c("#C84D4C", "#2883B1", "#ECB41F"),
               queries = list(
                 list(query = elements, params = list("Omic", "Metabolomic"), color = "#E1D2B2", active = TRUE, query.name = "Metabolomic"),
                 list(query = elements, params = list("Omic", "Proteomic"), color = "#A68A6D", active = TRUE, query.name = "Proteomic")
               ))

# Save PDF version of Figure 1
pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure1.pdf"),
    width = 180 / 25.4, height = 135 / 25.4, family = "Arial")
fig1b
dev.off()
