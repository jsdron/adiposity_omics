#########################################
# Script: tissue_express.R
# Description: ...
# Key Outputs:
#   - ...
#########################################

###### LIBRARIES ######
library(data.table)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(biomaRt)
library(UpSetR)

###### eQTL DATAFRAME ######
# Load baseline covariates 
main <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/target_protein.all_tissues_eqtls.tsv.gz")

main$fat_depot <- NA

main$fat_depot[main$gene_symbol == "THBS2"] <- "ASAT"

main$fat_depot[main$gene_symbol %in% c("ABL1", "CCL17", "MSR1", "CFB", "KHK", 
                                       "TSPAN8", "SELPLG", "ANXA2", "ASAH2", 
                                       "ADAMTSL5", "SHBG", "ITGB6")] <- "VAT"

main$fat_depot[main$gene_symbol %in% c("SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC")] <- 
  ifelse(is.na(main$fat_depot[main$gene_symbol %in% c("SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC")]), 
         "GFAT", 
         "VAT,GFAT")

### Convert rsIDs from MR IVs to GRCh38
# Connect to Ensembl SNP database (GRCh38)
# snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
# 
# # Read your rsID list
# rsids <- readLines("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/rsIDs.txt")  # One rsID per line
# 
# # Split into batches (e.g., 100 at a time)
# rsid_batches <- split(rsids, ceiling(seq_along(rsids) / 100))
# 
# # Query in batches
# results <- lapply(rsid_batches, function(batch) {
#   getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
#         filters = "snp_filter",
#         values = batch,
#         mart = snp_mart)
# })
# 
# # Combine and format
# snp_info <- do.call(rbind, results)
# snp_info <- snp_info %>%
#   dplyr::mutate(chr_pos = paste0("chr", chr_name, ":", chrom_start)) %>%
#   dplyr::rename(rsID = refsnp_id)
# 
# write.table(snp_info, "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/rsid_to_GRCh38_positions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

###### SIGNIFICANT eQTLS ######
threshold <- 0.05 / (18*50) # candidate proteins * tissue types

sig <- subset(main, main$pval_nominal < threshold)

### Restrict to the eQTLs used in the MR
rsids <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/rsid_to_GRCh38_positions.tsv",
               select = c("rsID", "chr_pos"))

sig$chr_pos <- paste0(sig$chr, ":", sig$pos)
sig <- subset(sig, sig$chr_pos %in% rsids$chr_pos)

table(sig$gene_symbol)
table(sig$tissue)
table(sig$gene_symbol, sig$tissue)

###### BAR PLOT######
plot_df <- sig %>%
  count(fat_depot, tissue, gene_symbol, fat_depot)

# Plot
ggplot(plot_df, aes(x = tissue, y = n, fill = gene_symbol)) +
  geom_bar(stat = "identity") +
  labs(x = "Tissue Type", y = "Count of Significant eQTLs", fill = "Gene Symbol") +
  theme_bw() +
  facet_wrap("fat_depot", nrow = 4, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot_df, aes(x = gene_symbol, y = n, fill = tissue)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(x = "Count of Significant eQTLs", y = "Tissue Type", fill = "Tissue Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot_df, aes(x = n, y = gene_symbol, fill = tissue)) +
  geom_bar(stat = "identity") +
  labs(x = "Count of Significant eQTLs", y = "Gene", fill = "Tissue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

upset_input <- main %>%
  filter(!is.na(gene_symbol), !is.na(tissue)) %>%
  distinct(gene_symbol, tissue) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = tissue, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_symbol")




mat <- table(sig$gene_symbol, sig$tissue)
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, 
         color = colorRampPalette(c("white", "red"))(100),
         display_numbers = TRUE, fontsize = 8)

