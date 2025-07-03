#########################################
# Script: pathway_figure.R
# Description: Generates pathway enrichment plots for proteins associated 
#              with ASAT, VAT, GFAT, and their overlaps.
# Key Outputs:
#   - Figure 3 (pathway enrichment plots)
#########################################

###### LIBRARIES ######
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(googlesheets4)
library(ggrepel)
library(gprofiler2)
library(ggrepel)
library(ggpubr)
library(scales)
library(extrafont)

font_import()  # This imports all available fonts to R
loadfonts(device = "pdf")  # Load fonts for PDF output

###### SET VARIABLES ######
# Define custom colors and shapes
custom_colors <- c("ASAT" = "#ECB41F", "VAT" = "#C84D4C", "GFAT" = "#2883B1") # Adjust colors as needed

# Date the analysis was completed 
input_date <- "20250127" 

# Correction factor
threshold = 41+445 

###### FUNCTIONS ######
# Lighten a color by blending with white
lighten_color <- function(hex, factor = 0.7) {
  rgb_color <- col2rgb(hex) / 255
  light_rgb <- rgb_color + (1 - rgb_color) * factor
  rgb(light_rgb[1], light_rgb[2], light_rgb[3])
}

# Run gProfiler for significant results only
gostres_sig_analysis <- function(i) {
  gost(query = query_list[i],
       organism = "hsapiens",
       domain_scope = "custom",
       custom_bg = olink_background,
       user_threshold = 0.05,
       correction_method = "fdr",
       significant = TRUE,
       sources = c("GO:BP", "KEGG", "REAC", "WP"),
       exclude_iea = FALSE,
       multi_query = FALSE)
}

# Run gProfiler for all results (used for full plot)
gostres_analysis <- function(i) {
  gost(query = query_list[i],
       organism = "hsapiens",
       domain_scope = "custom",
       custom_bg = olink_background,
       user_threshold = 0.05,
       correction_method = "fdr",
       significant = FALSE,
       sources = c("GO:BP", "KEGG", "REAC", "WP"),
       exclude_iea = FALSE,
       multi_query = FALSE)
}

# Build enrichment plot
gostres_plot <- function(gostres, label_data, colour_pick, plot_data) {
  lighter_colour <- lighten_color(colour_pick, factor = 0.2)

  plot <- gostplot(gostres, capped = FALSE, interactive = FALSE,
                   pal = c(`GO:BP` = "grey90", KEGG = "grey90", REAC = "grey90", WP = "grey90")) +
    labs(x = "Annotation", y = expression(-log[10]("P"))) +
    facet_wrap("query", scales = "free", nrow = 3) +
    scale_y_continuous(limits = c(0, 5)) +
    theme_cowplot() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 14, face = "bold", hjust = 0),
      strip.background = element_rect(fill = "white", colour = "white"),
      legend.text = element_text(size = 7),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      axis.ticks.x = element_blank()
    ) +
    geom_point(data = subset(plot_data, significant == TRUE),
               aes(x = order, y = logpval, size = term_size_scaled, alpha = 0.6),
               color = lighter_colour, show.legend = FALSE) +
    geom_point(data = label_data,
               aes(x = order, y = logpval, size = term_size_scaled, alpha = 1),
               color = colour_pick, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25, color = "black") +
    geom_text_repel(data = label_data,
                    aes(x = order, y = logpval, label = term_name),
                    size = 3, ylim = c(-log10(0.05), 5),
                    box.padding = 0.35,
                    point.padding = 0.3,
                    segment.color = 'black', segment.size = 0.3,
                    max.overlaps = Inf,
                    force = 2, force_pull = 1, nudge_y = 0.1)
  return(plot)
}

###### LOAD DATAFRAME(S) ######
# Load in association data and extract it down to the significant associations
output_all <- fread(paste0('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n5023.', input_date, '.tsv.gz'))
output_all$sig <- ifelse(output_all$p.value < 0.05 / threshold, "sig", "ns")
output_all$adi_label <- ifelse(output_all$outcome == "vatadjbmi3", "VAT",
                               ifelse(output_all$outcome == "asatadjbmi3", "ASAT", "GFAT"))

# The list of proteins covered by the Olink panel
olink_background <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/input/DAVID/olink_proteins.txt', header = TRUE)
olink_background <- olink_background$olink

###### LIST OF SIGNIFICANT PROTEINS ######
# Proteins associated with each fat depot
vat <- subset(output_all, adi_label == "VAT" & sig == "sig")$term
gfat <- subset(output_all, adi_label == "GFAT" & sig == "sig")$term
asat <- subset(output_all, adi_label == "ASAT" & sig == "sig")$term

# Proteins associated with a single fat depot
vat_only <- setdiff(vat, union(gfat, asat))
gfat_only <- setdiff(gfat, union(vat, asat))
asat_only <- setdiff(asat, union(gfat, vat))

# Proteins associated with all 3 fat depots
all_3 <- Reduce(intersect, list(vat, gfat, asat))

###### GENERATE SUPPLEMENTAL TABLE OF ENRICHED PROTEIN PATHWAYS ######
# Documentation: https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
gostres_sig <- gost(query = list("VAT (all)" = vat, "ASAT (all)" = asat, "GFAT (all)" = gfat,
                                 "VAT (only)" = vat_only, "ASAT (only)" = asat_only, "GFAT(only)" = gfat_only,
                                 "Overlap (all)" = all_3),
                    organism = "hsapiens", 
                    ordered_query = FALSE, 
                    multi_query = FALSE, 
                    exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = (0.05), correction_method = "fdr", significant = TRUE, highlight = FALSE,
                    domain_scope = "custom", custom_bg = olink_background, 
                    numeric_ns = "", 
                    sources = c("GO:BP", "KEGG", "REAC", "WP"), 
                    as_short_link = FALSE)

sig_results <- gostres_sig$result

# Remove root terms
remove <- c("WP:000000", "KEGG:00000")
sig_results <- subset(sig_results, !(term_id %in% remove))

# Export 
write.table(sig_results, 
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/enriched_protein_pathways.tsv.gz")), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

###### FIGURE 3 ######
query_list <- list("VAT (all)" = vat,
                   "ASAT (all)" = asat,
                   "GFAT (all)" = gfat,
                   "VAT (only)" = vat_only,
                   "ASAT (only)" = asat_only,
                   "GFAT (only)" = gfat_only,
                   "Overlap (all)" = all_3)

# vat_all (all protein associated with VAT)
  gostres_sig <- gostres_sig_analysis(1)
  gostres_sig$result <- subset(sig_results, !(term_id %in% remove))
  
  gostres <- gostres_analysis(1)
  gostres$result <- subset(gostres$result, !(term_id %in% remove))
  
  static_plot <- gostplot(gostres, capped = FALSE, interactive = FALSE)
  plot_data <- static_plot$data
  extra_hits <- subset(plot_data, plot_data$significant == TRUE & (plot_data$term_name %in% c("chemotaxis", "chemokine-mediated signaling pathway"))) # & plot_data$precision>=0.10
  top_hits <- subset(plot_data, plot_data$significant == TRUE & plot_data$precision >= 0.10) |>
    dplyr::group_by(query) |>
    dplyr::arrange(desc(logpval)) |>
    dplyr::slice_head(n = 5)

  label_data <- rbind(top_hits, extra_hits)
  
  vat_all <- gostres_plot(gostres, label_data, "#C84D4C", plot_data)

# asat_all (all protein associated with ASAT)
  gostres_sig <- gostres_sig_analysis(2)
  gostres_sig$result <- subset(sig_results, !(term_id %in% remove))
  
  gostres <- gostres_analysis(2)
  gostres$result <- subset(gostres$result, !(term_id %in% remove))
  
  static_plot <- gostplot(gostres, capped = FALSE, interactive = FALSE)
  plot_data <- static_plot$data
  extra_hits <- subset(plot_data, plot_data$significant==TRUE & (plot_data$term_name %in% c())) # & plot_data$precision>=0.10
  top_hits <- subset(plot_data, plot_data$significant==TRUE & plot_data$precision>=0.10) |>
    dplyr::group_by(query) |>
    dplyr::arrange(desc(logpval)) |>
    dplyr::slice_head(n = 5)
  
  asat_all <- gostres_plot(gostres, top_hits, "#ECB41F", plot_data)

# gfat_all (all protein associated with GFAT)
  gostres_sig <- gostres_sig_analysis(3)
  gostres_sig$result <- subset(sig_results, !(term_id %in% remove))
  
  gostres <- gostres_analysis(3)
  gostres$result <- subset(gostres$result, !(term_id %in% remove))
  
  static_plot <- gostplot(gostres, capped = FALSE, interactive = FALSE)
  plot_data <- static_plot$data
  extra_hits <- subset(plot_data, plot_data$significant==TRUE & (plot_data$term_name %in% c("lipid homeostasis", "regulation of plasma lipoprotein particle levels", "brown fat cell differentiation"))) # & plot_data$precision>=0.10
  top_hits <- subset(plot_data, plot_data$significant==TRUE & plot_data$precision>=0.10) |>
    dplyr::group_by(query) |>
    dplyr::arrange(desc(logpval)) |>
    dplyr::slice_head(n = 5)
  
  label_data <- rbind(top_hits, extra_hits)
  
  gfat_all <- gostres_plot(gostres, label_data, "#2883B1", plot_data)

## all_3 (all protein associated with all 3 fat depots)
  gostres_sig <- gostres_sig_analysis(7)
  gostres_sig$result <- subset(sig_results, !(term_id %in% remove))
  
  gostres <- gostres_analysis(7)
  gostres$result <- subset(gostres$result, !(term_id %in% remove))
  
  static_plot <- gostplot(gostres, capped = FALSE, interactive = FALSE)
  plot_data <- static_plot$data
  label_data <- subset(plot_data, plot_data$significant==TRUE & 
                         plot_data$precision>=0.10)
  
  all_fat <- gostres_plot(gostres, label_data, "black", plot_data)
  

main_figure <- ggarrange(asat_all, gfat_all, vat_all, all_fat,
                         labels = c("A", "B", "C", "D"),
                         ncol = 1, nrow = 4)

pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure5.pdf"), 
    width = 180/25.4, height = 220/25.4, family = "Arial") 
main_figure
dev.off()
