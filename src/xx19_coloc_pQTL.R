#########################################
# Script: coloc_pQTL.R
# Description: Colocalization analysis of pQTLs and fat depot GWAS
#   - Targets: MSR1–VAT, NFASC–GFAT
#   - Uses coloc.abf() with custom input preparation
#   - Generates posterior plots
#########################################

###### LIBRARIES ######
if(!require("remotes"))
  install.packages("remotes")

install_github("chr1swallace/coloc",build_vignettes= TRUE)

library(data.table)
library(coloc)
library(dplyr)
library(ggplot2)
library(stringr)

###### FUNCTIONS ######
### YOU'LL WANT TO DOUBLE CHECK THE COLUMN NAMES USED IN THIS FUNCTION! i know some of the inputs are column names, but still
prepare_coloc_inputs <- function(pqtl_path, gwas_path, gene_chr, region_start, region_end, 
                                 pqtl_snp_col = "rsid", gwas_snp_col = "rsid",
                                 pqtl_pos_col = "pos", gwas_pos_col = "pos") {
  
  # Read files
  pqtl <- fread(pqtl_path)
  gwas <- fread(gwas_path)
  
  # Filter to cis-region (± window)
  pqtl <- pqtl[pqtl[[pqtl_pos_col]] >= region_start & pqtl[[pqtl_pos_col]] <= region_end & pqtl$chr == gene_chr]
  gwas <- gwas[gwas[[gwas_pos_col]] >= region_start & gwas[[gwas_pos_col]] <= region_end & gwas$chr == gene_chr]
  
  # Create chr:pos SNP linker
  pqtl$snpid <- paste0(pqtl$chr, ":", pqtl[[pqtl_pos_col]])
  gwas$snpid <- paste0(gwas$chr, ":", gwas[[gwas_pos_col]])
  
  # Merge on SNPs in common
  common_snps <- intersect(pqtl$snpid, gwas$snpid)
  pqtl <- pqtl[pqtl$snpid %in% common_snps]
  gwas <- gwas[gwas$snpid %in% common_snps]
  
  # Harmonize effect direction if needed (simplified — assumes no strand issues)
  setkey(pqtl, snpid)
  setkey(gwas, snpid)
  merged <- merge(pqtl, gwas, by = "snpid", suffixes = c(".pqtl", ".gwas"))
  
  # Estimate varbeta if missing
  if (!"varbeta.pqtl" %in% names(merged)) {
    merged$varbeta.pqtl <- merged$SE.pqtl^2
  }
  if (!"varbeta.gwas" %in% names(merged)) {
    merged$varbeta.gwas <- merged$SE.gwas^2
  }
  
  # Add position column (needed for plot_dataset)
  merged$position <- merged[[paste0(pqtl_pos_col, ".pqtl")]]
  
  # Format datasets
  data1 <- list(
    snp = merged$snpid,
    beta = merged$beta.pqtl,
    varbeta = merged$varbeta.pqtl,
    MAF = merged$MAF.pqtl,
    N = unique(merged$N.pqtl),
    type = "quant",
    position = merged$position,
    pvalues = merged$pval.pqtl
  )
  
  data2 <- list(
    snp = merged$snpid,
    beta = merged$beta.gwas,
    varbeta = merged$varbeta.gwas,
    MAF = merged$MAF.gwas,
    N = unique(merged$N.gwas),
    type = "quant",
    position = merged$position,
    pvalues = merged$pval.gwas
  )
  
  return(list(data1 = data1, data2 = data2, merged = merged))
}

###### SET ANALYSIS PARAMETERS ######
# All in GRCh38
results_dir <- "results/coloc" # or wherever makes most sense
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

### THESE NEED TO BE FIXED FOR MSR1 AND NFASC! IM NO SURE WHAT THEY ARE! you can check from the MR script to see how it was done / where the chr, start, and end values came from
targets <- list(
  list(gene = "MSR1", depot = "VAT", chr = 8, start = 16000000, end = 16800000),
  list(gene = "NFASC", depot = "GFAT", chr = 1, start = 204000000, end = 205000000)
)

# File path templates - CAN GET FROM MR SCRIPT! 
pqtl_base <- "/path/to/pqtl/summarystats/"  
gwas_base <- "/path/to/fat_depot/gwas/"     

###### RUN coloc FOR EACH TARGET ######
for (t in targets) {
  gene <- t$gene
  depot <- t$depot
  chr <- t$chr
  region_start <- t$start
  region_end <- t$end
  
  pqtl_file <- paste0(pqtl_base, gene, "_cis_pqtl.txt.gz")     # TODO: verify filename
  gwas_file <- paste0(gwas_base, depot, "_GWAS_sumstat.txt.gz") # TODO: verify filename
  
  cat("Running coloc for", gene, "and", depot, "\n")
  
  coloc_inputs <- prepare_coloc_inputs(
    pqtl_path = pqtl_file,
    gwas_path = gwas_file,
    gene_chr = chr,
    region_start = region_start,
    region_end = region_end,
    pqtl_snp_col = "pos",   # adjust if needed
    gwas_snp_col = "pos"
  )
  
  d1 <- coloc_inputs$data1
  d2 <- coloc_inputs$data2
  merged <- coloc_inputs$merged
  
  check_dataset(d1)
  check_dataset(d2)
  
  res <- coloc.abf(dataset1 = d1, dataset2 = d2)
  summary_file <- paste0(results_dir, "/coloc_summary_", gene, "_", depot, ".txt")
  full_file <- paste0(results_dir, "/coloc_full_", gene, "_", depot, ".tsv")
  
  writeLines(capture.output(res$summary), summary_file)
  fwrite(res$results, full_file, sep = "\t")
  
  # Extract H4 SNPs
  h4_hits <- res$results[res$results$SNP.PP.H4 > 0.8, ]
  fwrite(h4_hits, paste0(results_dir, "/coloc_H4gt80_", gene, "_", depot, ".tsv"), sep = "\t")
  
  # Plot
  pdf(file = paste0(results_dir, "/coloc_plot_", gene, "_", depot, ".pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot_dataset(d1, main = paste(gene, "(pQTL)"))
  plot_dataset(d2, main = paste(depot, "(GWAS)"))
  dev.off()
}
