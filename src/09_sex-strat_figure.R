##############################################
# Script: sex-strat_figure.R
# Description: Generates volcano and sex-stratified scattter plots 
#              for protein and metabolite associations.
# Key Outputs:
#   - Supplemental Figure 5 (volcano and scatter plots)
##############################################

###### LIBRARIES ######
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(stats)

###### SET VARIABLES ######
# Correction factor
threshold <- 41+445 

# Nominal significance
nominal    <- -log10(0.05)

# Bonferroni corrected significance 
threshold2 <- -log10(0.05 / threshold)

###### PREPARE DATAFRAME(S) ######
# Load in protein association results for all sexes
prot_epi <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.20250127.tsv.gz")
prot_epi$Stratification.Group <- "All"

# Load in sex-stratified protein association results and merge
tmp <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.sex_stratified.n5023.20250127.tsv.gz")
prot_epi <- rbind(prot_epi, tmp)
prot_epi <- as.data.frame(prot_epi)

# Label fat depots and define plotting groups
prot_epi$outcome <- ifelse(prot_epi$outcome == "vatadjbmi3", "VAT",
                    ifelse(prot_epi$outcome == "asatadjbmi3", "ASAT", "GFAT"))
prot_epi$label <- paste0(prot_epi$term, "_", prot_epi$outcome)
prot_epi$alpha_value <- ifelse(prot_epi$p.value <= 0.05 / threshold, 1, 0.5)
prot_epi$sig <- ifelse(prot_epi$p.value <= 0.05 / threshold, "sig", 
                ifelse(prot_epi$p.value < 0.05 & prot_epi$p.value > 0.05 / threshold, "nominal", "ns"))
prot_epi$sig2 <- ifelse(prot_epi$p.value <= 0.05 / threshold, "sig", "ns")
prot_epi$group <- ifelse(prot_epi$beta < 0 & prot_epi$sig == "sig", "pos",
                  ifelse(prot_epi$beta > 0 & prot_epi$sig == "sig", "neg", "nothing"))

# Load in metabolomic association results for all sexes
metab_epi <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n22630.20250127.tsv.gz")
metab_epi$Stratification.Group <- "All"

# Load in sex-stratified metabolomic association results and merge
tmp <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.sex_stratified.n22630.20250127.tsv.gz")
metab_epi <- rbind(metab_epi, tmp)
metab_epi <- as.data.frame(metab_epi)

# Label fat depots and define plotting groups
metab_epi$label <- paste0(metab_epi$term, "_", metab_epi$outcome)
metab_epi$alpha_value <- ifelse(metab_epi$p.value <= 0.05 / threshold, 1, 0.5)
metab_epi$sig <- ifelse(metab_epi$p.value <= 0.05 / threshold, "sig", 
                 ifelse(metab_epi$p.value < 0.05 & metab_epi$p.value > 0.05 / threshold, "nominal", "ns"))
metab_epi$sig2 <- ifelse(metab_epi$p.value <= 0.05 / threshold, "sig", "ns")
metab_epi$group <- ifelse(metab_epi$beta < 0 & metab_epi$sig == "sig", "pos",
                   ifelse(metab_epi$beta > 0 & metab_epi$sig == "sig", "neg", "nothing"))

###### FIGURE PANELS A & B: VOLCANO PLOTS ######
# Proteomic volcano
b <- ggplot(prot_epi, aes(x = beta, y = -log10(p.value), color = group)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#dddddd") +
  geom_point(size = 1, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c("darkred", "lightgray", "darkblue")) +
  geom_hline(yintercept = threshold2, linewidth = 0.25, linetype = "dashed", color = "black") +
  labs(x = "Effect estimate", y = expression(-log[10]("P"))) +
  facet_grid(Stratification.Group ~ outcome, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 8),
    strip.text = element_text(size = 12)
  )

# Metabolomic volcano
a <- ggplot(metab_epi, aes(x = beta, y = -log10(p.value), color = group)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#dddddd") +
  geom_point(size = 1, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c("darkred", "lightgray", "darkblue")) +
  geom_hline(yintercept = threshold2, linewidth = 0.25, linetype = "dashed", color = "black") +
  labs(x = "Effect estimate", y = expression(-log[10]("P"))) +
  facet_grid(Stratification.Group ~ outcome, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )

sFig5ab <- ggarrange(a, b, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "bottom")

###### FIGURE PANEL C: METAB SCATTER ######
int_metab <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.sex_stratified.n22630.20250127.tsv.gz")
int_metab <- subset(int_metab, int_metab$sig_interaction == "Significant")
int_metab$linker <- paste0(int_metab$analyte, "_", int_metab$adi_label)
mr_metab <- c("HDL_CE", "M_HDL_CE", "M_HDL_C", "HDL_L", "HDL_PL", "M_HDL_FC", "M_HDL_P", "ApoA1", "M_HDL_L", "HDL_P", "M_HDL_PL")

wide_metab <- reshape(metab_epi[metab_epi$Stratification.Group %in% c("Female", "Male"), ],
                timevar = "Stratification.Group", idvar = c("term", "outcome"),
                direction = "wide")
wide_metab$group <- ifelse(wide_metab$sig2.Female=="sig" & wide_metab$sig2.Male=="sig", "Both",
                     ifelse(wide_metab$sig2.Female=="sig" & wide_metab$sig2.Male=="ns", "Female only",
                            ifelse(wide_metab$sig2.Female=="ns" & wide_metab$sig2.Male=="sig", "Male only", "Neither"))) 

wide_metab <- reshape(metab_epi[metab_epi$Stratification.Group %in% c("Female", "Male"), ],
                      timevar = "Stratification.Group", idvar = c("term", "outcome"),
                      direction = "wide")

wide_metab$group <- ifelse(wide_metab$sig2.Female == "sig" & wide_metab$sig2.Male == "sig", "Both",
                    ifelse(wide_metab$sig2.Female == "sig" & wide_metab$sig2.Male == "ns", "Female only",
                    ifelse(wide_metab$sig2.Female == "ns" & wide_metab$sig2.Male == "sig", "Male only", "Neither")))

metab_eff <- ggplot(wide_metab, aes(x = beta.Female, y = beta.Male, color = group)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#dddddd") +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#dddddd") +
  geom_point(size = 1, show.legend = TRUE) +
  geom_smooth(method = 'lm', formula = y ~ x, aes(group = 1), colour = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Neither" = "lightgray", "Both" = "#693182", 
                                "Male only" = "darkblue", "Female only" = "darkred")) +
  facet_grid(~outcome, scales = "free") +
  labs(x = "Female effect estimate", y = "Male effect estimate", color = "Significance") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 8),
    strip.text = element_text(size = 7)
  ) +
  ggpubr::stat_cor(r.accuracy = 0.01, cor.coef.name = "r", position = "identity",
                   aes(x = beta.Female, y = beta.Male, 
                       label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
                   inherit.aes = FALSE) +
  geom_text_repel(data = subset(wide_metab, wide_metab$group != "Neither" & 
                                label.Female %in% int_metab$linker & term %in% mr_metab), 
                  aes(x = beta.Female, y = beta.Male, label = term), 
                  size = 3, box.padding = 0.5, force = 10, 
                  segment.color = "black", segment.size = 0.3, 
                  show.legend = FALSE)

###### FIGURE PANEL D: PROT SCATTER ######
int_prot <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.sex_stratified.n5023.20250127.tsv.gz")
int_prot <- subset(int_prot, int_prot$sig_interaction == "Significant")
int_prot$linker <- paste0(int_prot$analyte, "_", int_prot$adi_label)

candidate_prot <- c("THBS2", "SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC", "ABL1", 
                    "CCL17", "MSR1", "CFB", "KHK", "TSPAN8", "SELPLG", "ANXA2", 
                    "ASAH2", "ADAMTSL5", "ITGB6")

wide_prot <- reshape(prot_epi[prot_epi$Stratification.Group %in% c("Female", "Male"), ],
                     timevar = "Stratification.Group", idvar = c("term", "outcome"),
                     direction = "wide")

wide_prot$group <- ifelse(wide_prot$sig2.Female == "sig" & wide_prot$sig2.Male == "sig", "Both",
                    ifelse(wide_prot$sig2.Female == "sig" & wide_prot$sig2.Male == "ns", "Female only",
                    ifelse(wide_prot$sig2.Female == "ns" & wide_prot$sig2.Male == "sig", "Male only", "Neither")))

prot_eff <- ggplot(wide_prot, aes(x = beta.Female, y = beta.Male, color = group)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#dddddd") +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#dddddd") +
  geom_point(size = 1, show.legend = TRUE) +
  geom_smooth(method = 'lm', formula = y ~ x, aes(group = 1), colour = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Neither" = "lightgray", "Both" = "#693182", 
                                "Male only" = "darkblue", "Female only" = "darkred")) +
  facet_grid(~outcome, scales = "free") +
  labs(x = "Female effect estimate", y = "Male effect estimate", color = "Significance") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 7)
  ) +
  ggpubr::stat_cor(r.accuracy = 0.01, cor.coef.name = "r", position = "identity",
                   aes(x = beta.Female, y = beta.Male, 
                       label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
                   inherit.aes = FALSE) +
  geom_text_repel(data = subset(wide_prot, wide_prot$group != "Neither" & 
                                label.Female %in% int_prot$linker & term %in% candidate_prot), 
                  aes(x = beta.Female, y = beta.Male, label = term), 
                  size = 3, box.padding = 0.5, force = 10, 
                  segment.color = "black", segment.size = 0.3, 
                  show.legend = FALSE)

sFig5cd <- ggarrange(metab_eff, prot_eff, labels = c("C", "D"), ncol = 1, common.legend = TRUE)

###### EXPORT PDF ######
pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/sFigure5.pdf", 
    width = 180 / 25.4, height = 10, family = "Arial")
ggarrange(sFig5ab, sFig5cd, nrow = 2, heights = c(2, 3))
dev.off()

###### TABLE OUTPUTS FOR MANUSCRIPT ######
table(wide_metab$group, wide_metab$outcome)
table(wide_prot$group, wide_prot$outcome)