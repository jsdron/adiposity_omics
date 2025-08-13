#########################################
# Script: MR_figure.R
# Description: Generates forward and reverse Mendelian Randomization (MR)
#              forest plots for metabolites and proteins associated with 
#              fat depots (ASAT, VAT, GFAT), including cross-directional effects.
# Key Outputs:
#   - Figure (composite MR forest plots for manuscript)
#########################################

###### LIBRARIES ######
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(grid)
library(tidytext)
library(colorspace)
library(extrafont)
loadfonts(device = "pdf")

###### FUNCTIONS ######
sort_and_collapse <- function(outcomes) {
  paste(sort(unique(outcomes)), collapse = ", ")
}

generate_palette <- function(base_color) {
  darker <- darken(base_color, amount = 0.5)
  lighter <- lighten(base_color, amount = 0.5)
  c(darker, base_color, lighter)
}

generate_palette2 <- function(base_color) {
  darker <- darken(base_color, amount = 0.5)
  lighter <- lighten(base_color, amount = 0.5)
  lightest <- lighten(lighter, amount = 0.5)
  c(darker, base_color, lighter, lightest)
}

###### SET VARIABLES ######
# Define colour palette
asat_colors <- generate_palette("#ECB41F")
gfat_colors <- generate_palette("#2883B1")
vat_colors  <- generate_palette("#C84D4C")

###### METABOLOMICS ######
# Load data
mr_all_metab <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/MR_2025-03-05/MR_metab_to_adipose.txt")
colnames(mr_all_metab) <- c("outcome","exposure","method","nsnp","b","se","pval","lo_ci","up_ci","or","or_lci95","or_uci95")

input_date <- "20250127"
epi_metab <- fread(paste0('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.', input_date, '.tsv.gz'),
                   select = c('term','estimate','std.error','statistic','p.value','outcome','adi_label'))

epi_metab$sig <- ifelse(epi_metab$p.value < 0.05 / (41 + 445), "sig", "ns")
epi_metab$adi_label <- ifelse(epi_metab$outcome == "vatadjbmi3", "VAT",
                              ifelse(epi_metab$outcome == "asatadjbmi3", "ASAT", "GFAT"))
epi_metab$label <- paste0(epi_metab$term, "_", epi_metab$adi_label)

# Preprocessing MR data
threshold <- 41
mr_all_metab$analysis <- ifelse(mr_all_metab$method %in% c("Inverse variance weighted (correlation inc)", "Wald ratio (correlation inc)",
                                                           "Inverse variance weighted", "Wald ratio"), "Primary", "Secondary")

mr_all_metab <- mr_all_metab %>%
  mutate(nominal_significance = ifelse(pval < 0.05, "Significant", "Not Significant"),
         significance = ifelse(pval < 0.05 / threshold, "Significant", "Not Significant"))

mr_all_metab$exposure <- gsub("met-d-", "", mr_all_metab$exposure)
mr_all_metab$label <- paste0(mr_all_metab$exposure, "_", mr_all_metab$outcome)
mr_all_metab <- subset(mr_all_metab, mr_all_metab$label %in% epi_metab$label[epi_metab$sig == "sig"])

# Forward MR
sig_IVW <- subset(mr_all_metab, analysis == "Primary" & nominal_significance == "Significant")
sig_df <- subset(mr_all_metab, label %in% sig_IVW$label)
sig_df <- subset(sig_df, !method %in% c("Egger intercept (correlation inc)", "Simple mode", "Weighted median", "Weighted mode"))

robust_df <- sig_df %>%
  group_by(label) %>%
  filter(all(b > 0) | all(b < 0)) %>%
  ungroup()

robust_df$method <- recode(robust_df$method,
                           "Inverse variance weighted" = "IVW",
                           "MR Egger" = "MR Egger",
                           "MR Raps" = "MR-RAPS",
                           "Wald ratio" = "Wald Ratio")

epi_robust <- subset(epi_metab, label %in% robust_df$label)[, c("label", "estimate", "p.value")]
merge_df <- merge(robust_df, epi_robust, by = "label")
merge_df <- merge_df[sign(merge_df$estimate) == sign(merge_df$b), ]
merge_df$label2 <- paste0(merge_df$label, "_", merge_df$method)
merge_df <- merge_df %>% distinct(label2, .keep_all = TRUE)
merge_df$method <- factor(merge_df$method, levels = c("IVW", setdiff(unique(merge_df$method), "IVW")))
merge_df <- subset(merge_df, !(label %in% merge_df$label[merge_df$pval >= 0.05]))

# Plot forward
desired_order <- c("M_HDL_PL", "HDL_P", "M_HDL_L", "ApoA1", "M_HDL_P", "M_HDL_FC", 
                   "HDL_PL", "HDL_L", "M_HDL_C", "M_HDL_CE", "HDL_CE")

merge_df$exposure <- factor(merge_df$exposure, levels = desired_order)

GFAT_fwd <- merge_df %>% filter(outcome == "GFAT") %>%
  ggplot(aes(x = reorder_within(exposure, b, within = outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method)) +
  scale_x_reordered() +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = FALSE) +
  scale_color_manual(values = c("IVW" = gfat_colors[1], "MR Egger" = gfat_colors[2], "MR-RAPS" = gfat_colors[3])) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method",
       title = "Exposure: Metabolite\nOutcome: Fat depot") +
  theme_bw() +
  facet_wrap(~outcome, scales = "free", ncol = length(unique(robust_df$exposure))) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Reverse MR
rev_mr <- fread("/Volumes/medpop_esp2/mpan/Projects/Adiposity/Adiposity_Omics/results/MR/ADIPOSE_MET_twosample/all_outcomes_mr_results_Adipose-MET_unique2025-03-16.csv")
rev_mr$label <- paste0(rev_mr$id.outcome, "_", rev_mr$id.exposure)

# Subset down to candidates
rev_mr <- subset(rev_mr, label %in% merge_df$label)

rev_mr$significance <- ifelse(rev_mr$pval <= 0.05, "Significant", "Not Significant")
rev_mr <- subset(rev_mr, method %in% c("Inverse variance weighted", "MR Egger", "MR Raps"))

rev_mr$method <- recode(rev_mr$method,
                        "Inverse variance weighted" = "IVW",
                        "MR Egger" = "MR Egger",
                        "MR Raps" = "MR-RAPS",
                        "Wald ratio" = "Wald Ratio")

rev_mr$method <- factor(rev_mr$method, levels = c("IVW", setdiff(unique(rev_mr$method), "IVW")))

GFAT_rev <- rev_mr %>% filter(id.exposure == "GFAT") %>%
  ggplot(aes(x = reorder_within(id.outcome, b, within = id.outcome),
             y = b, ymin = lo_ci, ymax = up_ci,
             color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = TRUE) +
  scale_color_manual(values = c("IVW" = gfat_colors[1], "MR Egger" = gfat_colors[2], "MR-RAPS" = gfat_colors[3])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance",
       title = "Exposure: Fat depot\nOutcome: Metabolite") +
  theme_bw() +
  facet_wrap(~id.exposure, scales = "free") +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Combine metab_plots
metab_plots <- ggarrange(GFAT_fwd, GFAT_rev, ncol = 2, labels = c("A", "B"),
                         legend = "bottom", common.legend = TRUE,
                         font.label = list(size = 7, face = "bold", color = "black"))

###### PROTEOMICS ######
# New colour palettes with an extra colour
asat_colors <- generate_palette2("#ECB41F")
gfat_colors <- generate_palette2("#2883B1")
vat_colors  <- generate_palette2("#C84D4C")

# Load data
prot_linker <- fread("/Volumes/medpop_esp2/mpan/Projects/Adiposity/Adiposity_Omics/data/olink_protein_map_3k.tsv")
mr_all_prot <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/MR_nonoverlap_fixed/estimates/for_analysis/proteins/all_merged_MR_results.2025-03-05.tsv")

mr_all_prot$analysis <- ifelse(mr_all_prot$method %in% c("Inverse variance weighted (correlation inc)", "Wald ratio (correlation inc)", "Wald ratio"), "Primary", "Secondary")
mr_all_prot <- mr_all_prot %>% mutate(significance = ifelse(pval < 0.05, "Significant", "Not Significant"))
mr_all_prot$label <- paste0(mr_all_prot$id.exposure, "_", mr_all_prot$id.outcome)

input_date <- "20250127"
epi <- fread(paste0('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.', input_date, '.tsv.gz'))
epi$sig <- ifelse(epi$p.value < 0.05 / (41 + 445), "sig", "ns")
epi$adi_label <- ifelse(epi$outcome == "vatadjbmi3", "VAT",
                        ifelse(epi$outcome == "asatadjbmi3", "ASAT", "GFAT"))
epi$label <- paste0(epi$term, "_", epi$adi_label)

epi_sig <- subset(epi, epi$sig == "sig")
mr_all_prot <- subset(mr_all_prot, mr_all_prot$label %in% epi_sig$label)

# Forward MR
sig_IVW <- subset(mr_all_prot, significance == "Significant" & analysis == "Primary")
sig_df <- subset(mr_all_prot, label %in% sig_IVW$label)
sig_df <- subset(sig_df, method != "Egger intercept (correlation inc)")

robust_df <- sig_df %>%
  group_by(label) %>%
  filter(all(b > 0) | all(b < 0)) %>%
  ungroup()

robust_df$method <- ifelse(robust_df$method == "Inverse variance weighted (correlation inc)", "IVW",
                           ifelse(robust_df$method == "Egger (correlation inc)", "MR Egger",
                                  ifelse(robust_df$method == "MR Raps", "MR-RAPS",
                                         ifelse(robust_df$method == "Wald ratio", "Wald Ratio", "Other"))))

epi_robust <- epi[, c("label", "estimate", "p.value")]
merge_df <- merge(robust_df, epi_robust, by = "label")
merge_df <- merge_df[sign(merge_df$estimate) == sign(merge_df$b), ]
merge_df$label2 <- paste0(merge_df$label, "_", merge_df$method)
merge_df <- merge_df %>% distinct(label2, .keep_all = TRUE)
merge_df$method <- factor(merge_df$method, levels = c("IVW", setdiff(unique(merge_df$method), "IVW")))
merge_df <- subset(merge_df, !(label %in% merge_df$label[merge_df$pval >= 0.05]))

# Legend
dummy <- merge_df
dummy[5, "significance"] <- "Not Significant"

legend_plot <- ggplot(dummy, aes(x = id.exposure, y = b, color = method, shape = significance, fill = significance)) +
  geom_point() + theme_bw() +
  theme(legend.position = "bottom", legend.justification = "center", legend.box = "vertical", legend.text = element_text(size = 11), legend.margin = margin()) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  scale_color_manual(values = c("IVW" = vat_colors[1], "MR Egger" = vat_colors[2], "MR-RAPS" = vat_colors[3], "Wald Ratio" = vat_colors[4])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

combined_legend <- get_legend(legend_plot)

# Forward MR plots for each fat depot
VAT <- merge_df %>% filter(id.outcome == "VAT") %>%
  ggplot(aes(x = reorder_within(id.exposure, b, within = id.outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = TRUE) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  scale_color_manual(values = c("IVW" = vat_colors[1], "MR Egger" = vat_colors[2], "MR-RAPS" = vat_colors[3], "Wald Ratio" = vat_colors[4])) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  theme_bw() +
  facet_wrap(~id.outcome, scales = "free", drop = TRUE, ncol = length(unique(robust_df$id.exposure))) +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
        #plot.margin = margin(0, 0.1, 0.1, -0.4, "cm"))

GFAT <- merge_df %>% filter(id.outcome == "GFAT") %>%
  ggplot(aes(x = reorder_within(id.exposure, b, within = id.outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = FALSE) +
  scale_color_manual(values = c("IVW" = gfat_colors[1], "MR Egger" = gfat_colors[2], "MR-RAPS" = gfat_colors[3], "Wald Ratio" = gfat_colors[4])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  theme_bw() +
  facet_wrap(~id.outcome, scales = "free", drop = TRUE, ncol = length(unique(robust_df$id.exposure))) +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
       # plot.margin = margin(0, 0.1, 0.1, -0.4, "cm"))

ASAT <- merge_df %>% filter(id.outcome == "ASAT") %>%
  ggplot(aes(x = reorder_within(id.exposure, b, within = id.outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = FALSE) +
  scale_color_manual(values = c("IVW" = asat_colors[1], "MR Egger" = asat_colors[2], "MR-RAPS" = asat_colors[3], "Wald Ratio" = asat_colors[4])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  theme_bw() +
  facet_wrap(~id.outcome, scales = "free", drop = TRUE, ncol = length(unique(robust_df$id.exposure))) +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
        #plot.margin = margin(0, 0.1, 0.3, -0.4, "cm"))

plots_fwd <- annotate_figure(ggarrange(VAT, ggarrange(ASAT, GFAT, nrow = 2, heights = c(1,2)), ncol = 2),
                             top = text_grob("Exposure: Protein\nOutcome: Fat depot", face = "bold", size = 7))
print(plots_fwd)

# Reverse MR
rev_mr <- fread("/Volumes/medpop_esp2/mpan/Projects/Adiposity/Adiposity_Omics/results/MR/ADIPOSE_PROT_twosample/all_outcomes_mr_results_2025-03-12.csv")
rev_mr$label <- paste0(rev_mr$id.outcome, "_", rev_mr$id.exposure)
rev_mr <- subset(rev_mr, label %in% merge_df$label)

rev_mr$significance <- ifelse(rev_mr$pval <= 0.05, "Significant", "Not Significant")
rev_mr <- subset(rev_mr, method != "Egger intercept (correlation inc)")

rev_mr$method <- ifelse(rev_mr$method == "Inverse variance weighted (correlation inc)", "IVW",
                        ifelse(rev_mr$method == "Egger (correlation inc)", "MR Egger",
                               ifelse(rev_mr$method == "MR Raps", "MR-RAPS",
                                      ifelse(rev_mr$method == "Wald ratio", "Wald Ratio", ""))))

rev_mr$method <- factor(rev_mr$method, levels = c("IVW", setdiff(unique(rev_mr$method), "IVW")))

# Reverse MR plots
VAT_rev <- rev_mr %>% filter(id.exposure == "VAT") %>%
  ggplot(aes(x = reorder_within(id.outcome, b, within = id.outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = TRUE) +
  scale_color_manual(values = c("IVW" = vat_colors[1], "MR Egger" = vat_colors[2], "MR-RAPS" = vat_colors[3], "Wald Ratio" = vat_colors[4])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  theme_bw() +
  facet_wrap(~id.exposure, scales = "free", drop = TRUE, ncol = length(unique(robust_df$id.outcome))) +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
       # plot.margin = margin(0, 0.1, 0.1, -0.4, "cm"))

GFAT_rev <- rev_mr %>% filter(id.exposure == "GFAT") %>%
  ggplot(aes(x = reorder_within(id.outcome, b, within = id.outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = FALSE) +
  scale_color_manual(values = c("IVW" = gfat_colors[1], "MR Egger" = gfat_colors[2], "MR-RAPS" = gfat_colors[3], "Wald Ratio" = gfat_colors[4])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  theme_bw() +
  facet_wrap(~id.exposure, scales = "free", drop = TRUE, ncol = length(unique(robust_df$id.outcome))) +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
        #plot.margin = margin(0, 0.1, 0.1, -0.4, "cm"))

ASAT_rev <- rev_mr %>% filter(id.exposure == "ASAT") %>%
  ggplot(aes(x = reorder_within(id.outcome, b, within = id.outcome), y = b, ymin = lo_ci, ymax = up_ci, color = method, shape = significance, fill = significance)) +
  geom_pointrange(position = position_dodge(width = 0.7), show.legend = FALSE) +
  scale_color_manual(values = c("IVW" = asat_colors[1], "MR Egger" = asat_colors[2], "MR-RAPS" = asat_colors[3], "Wald Ratio" = asat_colors[4])) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 21)) +
  scale_fill_manual(values = c("Significant" = NA, "Not Significant" = "white")) +
  labs(y = "Effect estimate (95% CI)", x = "", color = "MR Method", fill = "Significance", shape = "Significance") +
  theme_bw() +
  facet_wrap(~id.exposure, scales = "free", drop = TRUE, ncol = length(unique(robust_df$id.outcome))) +
  scale_x_reordered() +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = 7),   # facet header
    axis.title   = element_text(size = 7),
    axis.text    = element_text(size = 6),
    plot.title   = element_text(hjust = 0.5, size = 7, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
        #plot.margin = margin(0, 0.1, 0.3, -0.4, "cm"))

plots_rev <- annotate_figure(ggarrange(VAT_rev, ggarrange(ASAT_rev, GFAT_rev, nrow = 2, heights = c(1,2)), ncol = 2),
                             top = text_grob("Exposure: Fat depot\nOutcome: Protein", face = "bold", size = 7))

print(plots_rev)

# Combine plots
prot_plots <- plot_grid(
  ggarrange(plots_fwd, plots_rev, ncol = 2, labels = c("C", "D"),
            font.label = list(size = 7, face = "bold", color = "black")),
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(prot_plots)

###### FINAL FIGURES ######
# Combine metab_plots (top) and prot_plots (bottom) → save as Figure5
Figure4 <- plot_grid(
  metab_plots,
  prot_plots,
  ncol = 1,
  rel_heights = c(1,1)
)

# Display
print(Figure4)

# Save to PDF
pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/pdf/Figure4_MR.pdf",
    width = 180 / 25.4, 
    height = 180 / 25.4, 
    pointsize = 7,    
    family = "Arial")
Figure4
dev.off()
