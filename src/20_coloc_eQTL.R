#########################################
# Script: coloc_eQTL.R
# Description: Colocalization analysis of eQTLs and fat depot GWAS.
#
# Key Outputs:
#   - Generates coloc output result files
#########################################

###### LIBRARIES ######
if(!require("remotes"))
  install.packages("remotes")

install_github("chr1swallace/coloc",build_vignettes= TRUE)

library(data.table)
library(coloc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(stringr)
library(extrafont)
loadfonts(device = "pdf")

###### FUNCTIONS ######
strip_version <- function(x) sub("\\..*$", "", x)

parse_gtex_variant <- function(dt, id_col = "variant_id") {
  # "chr9_129713144_G_A_b38" -> chr,pos,ref,alt
  parts <- tstrsplit(dt[[id_col]], "_", fixed = TRUE)
  dt[, chr := parts[[1]]]
  dt[, pos := as.integer(parts[[2]])]
  dt[, ref := parts[[3]]]
  dt[, alt := parts[[4]]]
  dt
}

flip_alleles <- function(beta, ref1, alt1, ref2, alt2) {
  # If alleles are swapped (A1/A2 vs A2/A1), flip sign; else NA (no change)
  swap <- ref1 == alt2 & alt1 == ref2
  flip <- ifelse(swap, -1, 1)
  list(beta = beta * flip, swapped = swap)
}

is_ambiguous <- function(a1, a2) paste0(a1, a2) %in% c("AT","TA","CG","GC")

# Helper to get region per gene (from your BED)
gene_region <- function(g) {
  r <- bed[V4 == g][1]
  list(chr = r$V1, start = as.integer(r$V2), end = as.integer(r$V3), gene_id = r$V5, gene_symbol = r$V4)
}

gene_region2 <- function(g) {
  r <- bed[V4 == g][1]
  list(chr = r$V1, start = as.integer(r$V2), end = as.integer(r$V3),
       gene_id = r$V5, gene_symbol = r$V4)
}

tissue_N <- function(tiss) {
  idx <- match(tiss, N_eqtl_all$Tissue)
  if (is.na(idx)) NA_real_ else as.numeric(N_eqtl_all$N_RNASeq_and_Genotyped_samples[idx])
}

annotate_one <- function(row) {
  # row is a single-row data.table
  g      <- row$gene; tiss <- row$tissue; trait <- row$gwas
  gr     <- gene_region2(g)
  N_eqtl <- tissue_N(tiss)
  if (is.na(N_eqtl)) return(NULL)
  
  eqtl_df <- eqtl_main[tissue == tiss & gene_symbol == g]
  if (!nrow(eqtl_df)) return(NULL)
  
  gwas_dt <- gwas_list[[trait]]
  N_g     <- if (length(gwas_N) == 1L) gwas_N else gwas_N[[trait]]
  
  inp <- try(prepare_coloc_inputs_eqtl(
    eqtl = eqtl_df, gwas = gwas_dt,
    gene_id = gr$gene_id, gene_symbol = g, tissue = tiss,
    chr = gr$chr, region_start = gr$start, region_end = gr$end,
    N_eqtl = N_eqtl, N_gwas = N_g
  ), silent = TRUE)
  if (inherits(inp, "try-error")) return(NULL)
  
  m <- inp$merged
  if (!nrow(m)) return(NULL)
  
  # p-values & z-scores
  m[, pval_eqtl := 2*pnorm(-abs(beta_eqtl/se_eqtl))]
  m[, pval_gwas := 2*pnorm(-abs(beta_gwas/se_gwas))]
  m[, z_eqtl := beta_eqtl/se_eqtl]
  m[, z_gwas := beta_gwas/se_gwas]
  
  # leads
  lead_gwas <- m[which.min(pval_gwas)]
  lead_eqtl <- m[which.min(pval_eqtl)]
  
  # overlap & correlation
  same_lead <- identical(lead_gwas$snpid, lead_eqtl$snpid)
  # correlation on shared non-missing SNPs
  z_cor <- suppressWarnings(cor(m$z_eqtl, m$z_gwas, use="complete.obs"))
  
  data.table(
    gene = g, gene_id = gr$gene_id, tissue = tiss, gwas = trait,
    chr = gr$chr, region_start = gr$start, region_end = gr$end,
    nsnps = length(inp$data1$snp),
    PP0 = row$PP0, PP1 = row$PP1, PP2 = row$PP2, PP3 = row$PP3, PP4 = row$PP4,
    min_p_eqtl = min(m$pval_eqtl, na.rm=TRUE),
    min_p_gwas = min(m$pval_gwas, na.rm=TRUE),
    
    lead_snp_gwas = lead_gwas$snpid,
    lead_gwas_p   = lead_gwas$pval_gwas,
    lead_gwas_beta= lead_gwas$beta_gwas,
    
    lead_snp_eqtl = lead_eqtl$snpid,
    lead_eqtl_p   = lead_eqtl$pval_eqtl,
    lead_eqtl_beta= lead_eqtl$beta_eqtl,
    
    same_lead_snp = same_lead,
    zscore_cor    = z_cor
  )
}

prepare_coloc_inputs_eqtl <- function(
    eqtl, gwas,
    gene_id = NULL, gene_symbol = NULL, tissue = NULL,
    chr, region_start, region_end,
    # eqtl columns (GTEx allpairs-like)
    eqtl_gene_id_col = "gene_id",
    eqtl_gene_symbol_col = "gene_symbol",
    eqtl_tissue_col = "tissue",
    eqtl_variant_id_col = "variant_id",
    eqtl_beta_col = "slope",
    eqtl_se_col = "slope_se",
    eqtl_maf_col = "af",
    N_eqtl,                      # REQUIRED
    # gwas columns (summary stats)
    gwas_chr_col = "CHR",
    gwas_pos_col = "GENPOS",
    gwas_ref_col = "ALLELE0",    # non-effect (BOLT style)
    gwas_alt_col = "ALLELE1",    # effect allele (BOLT style)
    gwas_beta_col = "BETA",
    gwas_se_col = "SE",
    gwas_maf_col = "A1FREQ",
    N_gwas,                      # REQUIRED
    drop_ambiguous = TRUE
) {
  stopifnot(!missing(N_eqtl), !missing(N_gwas))
  
  # 1) Parse GTEx variant -> alleles
  parse_gtex_variant(eqtl, id_col = eqtl_variant_id_col)
  
  # 2) Optional gene/tissue filters
  if (!is.null(gene_id)) {
    eqtl[[eqtl_gene_id_col]] <- as.character(eqtl[[eqtl_gene_id_col]])
    eqtl <- eqtl[strip_version(get(eqtl_gene_id_col)) == strip_version(gene_id)]
  }
  if (!is.null(gene_symbol)) {
    eqtl <- eqtl[get(eqtl_gene_symbol_col) == gene_symbol]
  }
  if (!is.null(tissue)) {
    eqtl <- eqtl[get(eqtl_tissue_col) == tissue]
  }
  
  # 3) Region filter (both datasets)
  #    NOTE: do not mutate gwas$CHR globally; create a local chr string.
  chr_str <- chr
  eqtl <- eqtl[chr == chr_str & pos >= region_start & pos <= region_end]
  gwas <- gwas[get(gwas_chr_col) == sub("^chr","", chr_str) | get(gwas_chr_col) == chr_str]  # tolerate with/without "chr"
  gwas <- gwas[get(gwas_pos_col) >= region_start & get(gwas_pos_col) <= region_end]
  
  # 4) Build join key by position (robust to multi-allelic; alleles harmonized next)
  eqtl[, snp_pos := paste(chr, pos, sep=":")]
  gwas[, snp_pos := paste(
    ifelse(grepl("^chr", get(gwas_chr_col)), get(gwas_chr_col), paste0("chr", get(gwas_chr_col))),
    get(gwas_pos_col), sep=":"
  )]
  
  common_pos <- intersect(eqtl$snp_pos, gwas$snp_pos)
  if (!length(common_pos)) stop("No overlapping positions in region.")
  eqtl <- eqtl[snp_pos %in% common_pos]
  gwas <- gwas[snp_pos %in% common_pos]
  
  # 5) Merge on position; bring alleles and summary stats
  setkey(eqtl, snp_pos); setkey(gwas, snp_pos)
  m <- merge(
    eqtl[, .(snp_pos, chr, pos,
             eqtl_ref = ref, eqtl_alt = alt,
             beta_eqtl = get(eqtl_beta_col),
             se_eqtl   = get(eqtl_se_col),
             maf_eqtl  = get(eqtl_maf_col))],
    gwas[, .(snp_pos,
             gwas_ref = get(gwas_ref_col), gwas_alt = get(gwas_alt_col),
             beta_gwas = get(gwas_beta_col),
             se_gwas   = get(gwas_se_col),
             maf_gwas  = get(gwas_maf_col))],
    by = "snp_pos", all = FALSE
  )
  
  # Drop incomplete
  m <- m[is.finite(beta_eqtl) & is.finite(se_eqtl) & is.finite(beta_gwas) & is.finite(se_gwas)]
  
  # 6) Harmonize alleles to eQTL ALT as the effect allele
  match_ok   <- m$eqtl_ref == m$gwas_ref & m$eqtl_alt == m$gwas_alt
  match_swap <- m$eqtl_ref == m$gwas_alt & m$eqtl_alt == m$gwas_ref
  
  if (any(match_swap)) {
    idx <- which(match_swap)
    m$beta_gwas[idx] <- -m$beta_gwas[idx]
    tmp <- m$gwas_ref[idx]; m$gwas_ref[idx] <- m$gwas_alt[idx]; m$gwas_alt[idx] <- tmp
    if (!all(is.na(m$maf_gwas[idx]))) m$maf_gwas[idx] <- 1 - m$maf_gwas[idx]
  }
  
  # Optionally drop ambiguous strand SNPs if alleles don’t align
  if (drop_ambiguous) {
    ambiguous <- is_ambiguous(m$eqtl_ref, m$eqtl_alt) | is_ambiguous(m$gwas_ref, m$gwas_alt)
  } else ambiguous <- rep(FALSE, nrow(m))
  
  keep <- (match_ok | match_swap) & !ambiguous
  m <- m[keep]
  
  # 7) Final unique key AFTER harmonization: chr:pos:ref:alt (using eQTL alleles)
  m[, snpid := paste(chr, pos, eqtl_ref, eqtl_alt, sep=":")]
  
  # 8) Variances and MAF
  m[, varbeta_eqtl := se_eqtl^2]
  m[, varbeta_gwas := se_gwas^2]
  m[, MAF := fifelse(is.finite(maf_gwas), maf_gwas,
                     fifelse(is.finite(maf_eqtl), maf_eqtl, NA_real_))]
  
  # 9) Deduplicate exact variants (keep smallest variance)
  setorder(m, snpid, varbeta_eqtl, varbeta_gwas)
  m <- m[, .SD[1], by = snpid]
  
  # 10) Build coloc datasets (ordered)
  setorder(m, chr, pos)
  data1 <- list( # eQTL
    snp      = m$snpid,
    beta     = m$beta_eqtl,
    varbeta  = m$varbeta_eqtl,
    MAF      = m$MAF,
    N        = N_eqtl,
    type     = "quant",
    position = m$pos
  )
  data2 <- list( # GWAS
    snp      = m$snpid,
    beta     = m$beta_gwas,
    varbeta  = m$varbeta_gwas,
    MAF      = m$MAF,
    N        = N_gwas,
    type     = "quant",
    position = m$pos
  )
  
  list(data1 = data1, data2 = data2, merged = m[])
}

###### RUN COLOC ######
# Read-in input files
eqtl_main <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/target_protein.all_tissues_eqtls.tsv.gz")
bed <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/gene_windows.bed")
N_eqtl_all <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/tissue_sample_size.txt")

gwas_paths <- list(
  GFAT = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_gfatadjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz",
  ASAT = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_asatadjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz",
  VAT  = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_vatadjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz"
)

read_gwas <- function(p) {
  dt <- fread(p)
  # Ensure expected columns exist; rename if needed
  need <- c("CHR","GENPOS","ALLELE0","ALLELE1","BETA","SE","A1FREQ")
  has  <- toupper(names(dt))
  setnames(dt, names(dt), has)
  # Common alternates:
  altmap <- c("CHROM"="CHR","POS"="GENPOS","A0"="ALLELE0","A1"="ALLELE1",
              "BETA_SNP"="BETA","SE_SNP"="SE","MAF"="A1FREQ","AF"="A1FREQ")
  for (k in names(altmap)) if (k %in% names(dt) && !(altmap[k] %in% names(dt))) setnames(dt, k, altmap[k])
  stopifnot(all(need %in% names(dt)))
  dt[]
}

gwas_list <- lapply(gwas_paths, read_gwas)

# Sample sizes
gwas_N <- 32950

## Genes and tissues to scan
genes <- unique(bed$V4)                              
tissues <- unique(N_eqtl_all$Tissue)

## Main coloc loop (ABF only)
extract_abf_summary <- function(abf) {
  s <- abf$summary
  if (is.null(s)) return(NULL)
  if (is.null(dim(s))) {
    # named vector
    list(
      nsnps = as.integer(s[["nsnps"]]),
      PP0   = as.numeric(s[["PP.H0.abf"]]),
      PP1   = as.numeric(s[["PP.H1.abf"]]),
      PP2   = as.numeric(s[["PP.H2.abf"]]),
      PP3   = as.numeric(s[["PP.H3.abf"]]),
      PP4   = as.numeric(s[["PP.H4.abf"]])
    )
  } else {
    # 1-row data.frame/matrix
    list(
      nsnps = as.integer(s[1, "nsnps"]),
      PP0   = as.numeric(s[1, "PP.H0.abf"]),
      PP1   = as.numeric(s[1, "PP.H1.abf"]),
      PP2   = as.numeric(s[1, "PP.H2.abf"]),
      PP3   = as.numeric(s[1, "PP.H3.abf"]),
      PP4   = as.numeric(s[1, "PP.H4.abf"])
    )
  }
}

ix <- 1L
results <- list()
inp_dir <- "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/coloc"
if (!dir.exists(inp_dir)) dir.create(inp_dir, recursive = TRUE)

safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

for (g in genes) {
  gr <- gene_region(g)
  for (tiss in tissues) {
    
    # robust tissue-N lookup
    idx <- match(tiss, N_eqtl_all$Tissue)
    if (is.na(idx)) next
    N_eqtl <- as.numeric(N_eqtl_all$N_RNASeq_and_Genotyped_samples[idx])
    if (length(N_eqtl) == 0L || is.na(N_eqtl)) next
    
    # pre-slice eQTL; skip if empty
    eqtl_df <- subset(eqtl_main, tissue == tiss & gene_symbol == g)
    if (!nrow(eqtl_df)) next
    
    for (gw_name in names(gwas_list)) {
      gwas_dt <- gwas_list[[gw_name]]
      N_gwas  <- gwas_N
      if (length(N_gwas)==0L || is.na(N_gwas)) next
      
      inp <- try(prepare_coloc_inputs_eqtl(
        eqtl = eqtl_df, gwas = gwas_dt,
        gene_id = gr$gene_id, gene_symbol = gr$gene_symbol, tissue = tiss,
        chr = gr$chr, region_start = gr$start, region_end = gr$end,
        N_eqtl = N_eqtl, N_gwas = N_gwas
      ), silent = TRUE)
      if (inherits(inp, "try-error")) next
      
      n_snps <- length(inp$data1$snp)
      if (n_snps < 50L) next
      
      abf <- try(coloc.abf(inp$data1, inp$data2), silent = TRUE)
      if (inherits(abf, "try-error")) next
      
      s <- extract_abf_summary(abf)   # your helper from earlier
      if (is.null(s)) next
      
      m <- inp$merged
      m[, pval_eqtl := 2*pnorm(-abs(beta_eqtl/se_eqtl))]
      m[, pval_gwas := 2*pnorm(-abs(beta_gwas/se_gwas))]
      
      fname <- sprintf("%s__%s__%s__%s_%d_%d.rds",
                       safe(gr$gene_symbol), safe(tiss), safe(gw_name),
                       safe(gr$chr), gr$start, gr$end)
      fpath <- file.path(inp_dir, fname)
      
      # store exactly what you’ll want to inspect later
      saveRDS(
        list(
          meta   = list(gene=gr$gene_symbol, gene_id=gr$gene_id, tissue=tiss, gwas=gw_name,
                        chr=gr$chr, start=gr$start, end=gr$end, N_eqtl=N_eqtl, N_gwas=N_gwas),
          data1  = inp$data1,     # eQTL input to coloc.abf
          data2  = inp$data2,     # GWAS input to coloc.abf
          merged = m              # your merged table with computed p-values
        ),
        fpath, compress = "xz"
      )
      
      results[[ix]] <- data.table(
        gene = gr$gene_symbol,
        gene_id = gr$gene_id,
        tissue = tiss,
        gwas = gw_name,
        chr = gr$chr, region_start = gr$start, region_end = gr$end,
        nsnps = s$nsnps,
        PP0 = s$PP0, PP1 = s$PP1, PP2 = s$PP2, PP3 = s$PP3, PP4 = s$PP4,
        min_p_eqtl = min(m$pval_eqtl, na.rm=TRUE),
        min_p_gwas = min(m$pval_gwas, na.rm=TRUE),
        inp_file = fpath
      )
      ix <- ix + 1L
    }
  }
}

# Save the results
res <- data.table::rbindlist(results, fill = TRUE)

# Define keep lists
keep_pairs <- list(
  ASAT = c("THBS2"),
  VAT  = c("ABL1", "CCL17", "MSR1", "CFB", "KHK", "TSPAN8", "SELPLG",
           "ANXA2", "ASAH2", "ADAMTSL5", "SHBG", "ITGB6"),
  GFAT = c("SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC")
)

# Filter res
res_filtered <- res[
  mapply(function(gwas_val, gene_val) {
    gwas_val %in% names(keep_pairs) &&
      gene_val %in% keep_pairs[[gwas_val]]
  }, gwas, gene)
]

data.table::fwrite(res, "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/tissue_coloc/coloc_abf_eQTL-GWAS.2025-08-26.tsv", sep = "\t")
data.table::fwrite(res_filtered, "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/tissue_coloc/coloc_abf_eQTL-GWAS_candidates.2025-08-26.tsv", sep = "\t")

# Save this file of regions for the LD matrix
loci_manifest <- res_filtered[, c(1,5:7)]
loci_manifest <- distinct(loci_manifest)
loci_manifest$chr <- gsub("chr", "", loci_manifest$chr)
colnames(loci_manifest) <- c("locus_id", "chr", "start_bp", "end_bp")

# Can stick with eQTL windows since they contain the pQTL windows
data.table::fwrite(loci_manifest, "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/locuszoom/loci_manifest.2025-08-26.tsv", sep = "\t")

###### VISUALIZE AND SUMMARIZE RESULTS ######
# inp_dir <- "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/coloc"
# 
# get_inp <- function(g, t, gw, dir = inp_dir) {
#   f <- file.path(dir, sprintf("%s_%s_%s.rds", g, t, gw))
#   if (!file.exists(f)) stop("Not found: ", f)
#   readRDS(f)
# }
# 
# gene <- "SHBG"
# tissue <- "Brain_Nucleus_accumbens_basal_ganglia"
# depot <- "GFAT"
# 
# inp <- get_inp(gene, tissue, depot)
# d1 <- inp$data1
# d2 <- inp$data2
# 
# plot_dataset(d1, main = paste0("eQTL (", tissue, ")"))
# plot_dataset(d2, main = paste0(depot, " (GWAS)"))







# # Overall PP4 distribution
# summary(res_filtered$PP4)
# hist(res_filtered$PP4, breaks=50, main="PP4 distribution", xlab="PP4")
# 
# res_filtered_modhigh <- res_filtered %>%
#   filter(PP3 > 0.6 | PP4 > 0.6)
# 
# # Melt PP3 and PP4 into long format for plotting
# df_eqtl_long <- res_filtered_modhigh %>%
#   select(gene, tissue, PP3, PP4, gwas) %>%
#   pivot_longer(cols = c(PP3, PP4), names_to = "PP_type", values_to = "PP_value") %>%
#   mutate(
#     PP_label = if_else(PP_value < 0.01, "<0.01", sprintf("%.2f", round(PP_value, 2))),
#     PP_type  = factor(PP_type, levels = c("PP3", "PP4"))
#   )
# 
# # Choose the tissue order you want (example puts adipose first, then others alphabetically)
# tissue_levels <- c(
#   "Adipose_Subcutaneous",
#   "Adipose_Visceral_Omentum",
#   sort(setdiff(unique(df_eqtl_long$tissue),
#                c("Adipose_Subcutaneous","Adipose_Visceral_Omentum")))
# )
# 
# # Build y-levels = (tissue1·PP3, tissue1·PP4, tissue2·PP3, tissue2·PP4, ...)
# pp_levels <- c("PP3","PP4")
# y_levels  <- as.vector(t(outer(tissue_levels, pp_levels, paste, sep = " · ")))
# 
# # Colour by tissue types -----
# df_eqtl_long2 <- df_eqtl_long %>%
#   mutate(
#     tissue  = factor(tissue, levels = tissue_levels),
#     PP_type = factor(PP_type, levels = pp_levels),
#     y_combo = factor(paste(tissue, PP_type, sep = " · "), levels = y_levels),
#     PP_label = ifelse(PP_value < 0.01, "<0.01", sprintf("%.2f", round(PP_value, 2)))
#   )
# 
# # Colours that align with the GTEx tissues 
# tissue_pal <- c(
#   "Adipose_Subcutaneous"                 = "#FF6601",
#   "Artery_Tibial"                        = "#FF0001",
#   "Brain_Nucleus_accumbens_basal_ganglia"= "#EEEE00",
#   "Lung"                                 = "#99FE01",
#   "Adipose_Visceral_Omentum"             = "#FFAA00",
#   "Adrenal_Gland"                        = "#32DD33",
#   "Brain_Cortex"                         = "#EEEE00",
#   "Cells_EBV-transformed_lymphocytes"    = "#CC66FF",
#   "Esophagus_Mucosa"                     = "#552200",
#   "Pituitary"                            = "#AAFF99",
#   "Testis"                               = "#AAAAAA"
# )
# 
# ggplot(
#   df_eqtl_long2,
#   aes(x = PP_value, y = y_combo, fill = tissue)
# ) +
#   annotate("rect", xmin = 0.6, xmax = 0.8, ymin = -Inf, ymax = Inf,
#            fill = "khaki", alpha = 0.25, color = NA) +
#   annotate("rect", xmin = 0.8, xmax = 1.0, ymin = -Inf, ymax = Inf,
#            fill = "lightgreen", alpha = 0.25, color = NA) +
#   geom_col(width = 0.7) +
#   geom_text(aes(x = pmin(PP_value + 0.01, 1.0), label = PP_label),
#             hjust = 0, size = 2.2, color = "black") +
#   facet_grid(rows = vars(gene), scales = "free_y", space = "free_y") +
#   scale_y_discrete(expand = expansion(add = c(0.6, 0.6), mult = c(0, 0))) +
#   scale_x_continuous(limits = c(0, 1.03),
#                      expand = expansion(mult = 0, add = c(0.01, 0))) +
#   labs(x = "Posterior Probability", y = "", fill = "Tissue",
#        title = "PP3 vs. PP4 from eQTL colocalization analysis") +
#   scale_y_discrete(labels = function(x) sub("^.*·\\s*", "", x)) +
#   theme_bw(base_size = 7) +
#   theme(
#     axis.text.x  = element_text(size = 6),
#     axis.text.y  = element_text(size = 6),
#     axis.title.x = element_text(size = 7),
#     axis.title.y = element_blank(),
#     strip.text   = element_text(size = 7, margin = margin(1,1,1,1)),
#     legend.text  = element_text(size = 7),
#     legend.title = element_text(size = 7),
#     plot.title = element_text(size = 7, face = "bold"),
#     legend.position   = "bottom",
#     legend.key.height = unit(0.3, "lines"),
#     legend.key.width  = unit(0.6, "lines"),
#     legend.spacing.x  = unit(0.1, "cm"),
#     legend.spacing.y  = unit(0.1, "cm"),
#     panel.grid.major.y = element_blank(),
#     panel.spacing.y     = unit(3, "pt"),
#     panel.spacing.x     = unit(0, "pt"),
#     plot.margin         = margin(2, 2, 2, 2)
#   )
# 
# # Colour by PP_types -----
# df_eqtl_long2 <- df_eqtl_long %>%
#   mutate(
#     tissue  = factor(tissue, levels = tissue_levels),
#     PP_type = factor(PP_type, levels = pp_levels),
#     y_combo = factor(paste(tissue, PP_type, sep = " · "), levels = y_levels),
#     PP_label = ifelse(PP_value < 0.01, "<0.01", sprintf("%.2f", round(PP_value, 2)))
#   )
# 
# ggplot(
#   df_eqtl_long2,
#   aes(x = PP_value, y = y_combo, fill = PP_type)   # color by PP_type (as in pQTL)
# ) +
#   # Evidence bands
#   annotate("rect", xmin = 0.6, xmax = 0.8, ymin = -Inf, ymax = Inf,
#            fill = "khaki", alpha = 0.25, color = NA) +
#   annotate("rect", xmin = 0.8, xmax = 1.0, ymin = -Inf, ymax = Inf,
#            fill = "lightgreen", alpha = 0.25, color = NA) +
#   
#   # Bars + value labels
#   geom_col(width = 0.7) +
#   geom_text(aes(x = pmin(PP_value + 0.01, 1.0), label = PP_label),
#             hjust = 0, size = 2.2, color = "black") +
#   
#   # One facet per gene; consistent bar thickness
#   facet_grid(rows = vars(gene), scales = "free_y", space = "free_y") +
#   
#   # Show tissue names on the y-axis (remove " · PP3/PP4" suffix)
#   scale_y_discrete(
#     labels = function(x) sub("\\s*·\\s*PP[34]$", "", x),
#     expand = expansion(add = c(0.6, 0.6), mult = c(0, 0))
#   ) +
#   
#   # Tight x scale
#   scale_x_continuous(limits = c(0, 1.03),
#                      expand = expansion(mult = 0, add = c(0.01, 0))) +
#   
#   # PP_type colors (same as pQTL)
#   scale_fill_manual(
#     values = c("PP3" = "lightgrey", "PP4" = "#4AA1A1"),
#     breaks = c("PP3", "PP4"),
#     name   = "Hypothesis"
#   ) +
#   
#   labs(
#     x = "Posterior Probability",
#     y = "",
#     title = "PP3 vs. PP4 from eQTL colocalization analysis"
#   ) +
#   theme_bw(base_size = 7) +
#   theme(
#     axis.text.x  = element_text(size = 6),
#     axis.text.y  = element_text(size = 6),
#     axis.title.x = element_text(size = 7),
#     axis.title.y = element_blank(),
#     strip.text   = element_text(size = 7, margin = margin(1,1,1,1)),
#     legend.text  = element_text(size = 7),
#     legend.title = element_text(size = 7),
#     legend.position   = "bottom",
#     plot.title = element_text(size = 7, face = "bold"),
#     legend.key.height = unit(0.3, "lines"),
#     legend.key.width  = unit(0.6, "lines"),
#     legend.spacing.x  = unit(0.1, "cm"),
#     legend.spacing.y  = unit(0.1, "cm"),
#     panel.grid.major.y = element_blank(),
#     panel.spacing.y     = unit(3, "pt"),
#     panel.spacing.x     = unit(0, "pt"),
#     plot.margin         = margin(2, 2, 2, 2)
#   )

