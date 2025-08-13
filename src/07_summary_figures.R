#########################################
# Script: summary_figures.R
# Description: Combines metabolomic and proteomic results to:
#              (1) compare significant analytes across fat depots, 
#              (2) visualize shared vs. distinct signals via barplots and upset plots.
# Key Outputs:
#   - Figure (upset plot and barplots)
#########################################

###### LIBRARIES ######
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)   
library(ggpubr)    
library(UpSetR)
library(gridExtra)
library(extrafont)
library(ggplotify)

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

###### FIGURE: BARPLOTS OF SIGNIFICANT ANALYTES ######
plot_data <- output_all %>%
  count(adi_label, omic, sig) %>%
  group_by(adi_label, omic) %>%
  mutate(percentage = n / sum(n) * 100,
         rounded_per = round(percentage, 1))

# Absolute count plot
fig1a <- ggplot(plot_data[plot_data$sig == "sig", ], aes(x = adi_label, y = n, fill = omic)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(x = "Fat depot", y = "Significant analytes (n)", fill = "") +
  scale_fill_manual(values = c("#E1D2B2", "#A68A6D")) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.4)

# Percent bar plot
fig1b <- ggplot(plot_data[plot_data$sig == "sig", ], aes(x = adi_label, y = rounded_per, fill = omic)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  labs(x = "Fat depot", y = "Significant analytes (%)", fill = "") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  scale_fill_manual(values = c("#E1D2B2", "#A68A6D")) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)) +
  geom_text(aes(label = scales::percent(rounded_per / 100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 2.4)

# Combine for Supplemental Figure 2
Fig1ab <- ggarrange(fig1a, fig1b, ncol = 1, align = "h", labels = c("A", "B"),
                    font.label = list(size = 7, face = "bold", color = "black"))

svg(
  paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure1ab.svg"),
  width  = 2.5,  
  height = 4,
  family = "Arial"
)
Fig1ab
dev.off()

###### FIGURE: UPSET PLOT ######
sig_analyte <- output_all[output_all$sig == "sig", ]

# Get significant analytes by depot
A <- sig_analyte[sig_analyte$adi_label == "VAT", ]$term
B <- sig_analyte[sig_analyte$adi_label == "ASAT", ]$term
C <- sig_analyte[sig_analyte$adi_label == "GFAT", ]$term

# Build upset_df with one row per analyte and omic
upset_df <- sig_analyte %>%
  dplyr::distinct(Analyte = term, Omic = omic) %>%
  dplyr::mutate(
    VAT  = as.integer(Analyte %in% A),
    ASAT = as.integer(Analyte %in% B),
    GFAT = as.integer(Analyte %in% C)
  )

upset_df$VAT <- ifelse(upset_df$Analyte %in% A, 1, 0)
upset_df$ASAT <- ifelse(upset_df$Analyte %in% B, 1, 0)
upset_df$GFAT <- ifelse(upset_df$Analyte %in% C, 1, 0)

# Create upset plot
Fig1c <-  upset(upset_df, 
                sets = c("VAT", "ASAT", "GFAT"), 
                order.by = "freq",
                mainbar.y.label = "Intersection size", 
                sets.x.label = "Signficiant analytes",
                empty.intersections = "on",
                query.legend = "bottom",
                text.scale = rep(0.7, 6),  # Adjusts all text elements to 7pt.
                point.size = 3,
                sets.bar.color = c("#C84D4C", "#2883B1", "#ECB41F"),  # Change these colors as needed
                queries = list(list(query = elements, 
                                    params = list("Omic", c("Metabolomic", "Proteomic")), 
                                    color = "#E1D2B2", active = T, query.name = "Metabolomic"), #ababab
                               list(query = elements, 
                                    params = list("Omic", "Proteomic"), 
                                    color = "#A68A6D", active = T, query.name = "Proteomic")) #3f3f3f
)

pdf(
  paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure1c.pdf"),
  width  = 4.6,  
  height = 4,
  family = "Arial"
)
Fig1c
dev.off()
