#########################################
# Script: coloc_pQTL.R
# Description: Colocalization analysis of pQTLs and fat depot GWAS
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
norm_chr <- function(x){ x <- as.character(x); ifelse(grepl("^chr", x), x, paste0("chr", x)) }

is_ambiguous <- function(a1,a2) paste0(a1,a2) %in% c("AT","TA","CG","GC")

get_pqtl_N <- function(pqtl_file) {
  # read just the N column (fast); tolerate 'N' or 'n'
  hdr <- names(fread(pqtl_file, nrows = 0))
  Ncol <- if ("N" %in% hdr) "N" else if ("n" %in% hdr) "n" else stop("No 'N' column in pQTL file: ", pqtl_file)
  Nvec <- fread(pqtl_file, select = Ncol)[[1]]
  u <- unique(na.omit(Nvec))
  if (!length(u)) stop("All N are NA in pQTL file: ", pqtl_file)
  if (length(u) > 1) {
    warning("Multiple N values in ", pqtl_file, " (using median).")
    as.integer(round(median(u)))
  } else {
    as.integer(u[1])
  }
}

prepare_coloc_inputs_pqtl_file <- function(
    pqtl_path, gwas,
    chr, region_start, region_end,
    # pQTL columns
    pqtl_chr_col = "CHROM",
    pqtl_pos_col = "GENPOS",
    pqtl_ref_col = "ALLELE0",      # non-effect allele in pQTL
    pqtl_alt_col = "ALLELE1",      # effect allele in pQTL
    pqtl_beta_col= "BETA",
    pqtl_se_col  = "SE",
    pqtl_maf_col = "A1FREQ",
    N_pqtl,                         
    # GWAS columns 
    gwas_chr_col = "CHR",
    gwas_pos_col = "GENPOS",
    gwas_ref_col = "ALLELE0",
    gwas_alt_col = "ALLELE1",
    gwas_beta_col= "BETA",
    gwas_se_col  = "SE",
    gwas_maf_col = "A1FREQ",
    N_gwas,                        
    drop_ambiguous = TRUE
){
  stopifnot(!missing(N_pqtl), !missing(N_gwas))
  pqtl <- fread(pqtl_path)
  
  # 1) Region filter with tolerant chr formatting
  pqtl[, CHR_ := norm_chr(get(pqtl_chr_col))]
  gwas[, CHR_ := norm_chr(get(gwas_chr_col))]  # view; does not overwrite original columns
  chr_ <- norm_chr(chr)
  
  pqtl <- pqtl[CHR_ == chr_ & get(pqtl_pos_col) >= region_start & get(pqtl_pos_col) <= region_end]
  gwas <- gwas[CHR_ == chr_ & get(gwas_pos_col) >= region_start & get(gwas_pos_col) <= region_end]
  if (!nrow(pqtl) || !nrow(gwas)) stop("No rows after region filter.")
  
  # 2) Join by position first
  pqtl[, snp_pos := paste(CHR_, get(pqtl_pos_col), sep=":")]
  gwas[, snp_pos := paste(CHR_, get(gwas_pos_col), sep=":")]
  common_pos <- intersect(pqtl$snp_pos, gwas$snp_pos)
  if (!length(common_pos)) stop("No overlapping positions in region.")
  
  setkey(pqtl, snp_pos); setkey(gwas, snp_pos)
  m <- merge(
    pqtl[, .(snp_pos,
             chr = CHR_, pos = get(pqtl_pos_col),
             p_ref = get(pqtl_ref_col), p_alt = get(pqtl_alt_col),
             beta_p = get(pqtl_beta_col), se_p = get(pqtl_se_col),
             maf_p  = get(pqtl_maf_col))],
    gwas[, .(snp_pos,
             g_ref = get(gwas_ref_col), g_alt = get(gwas_alt_col),
             beta_g = get(gwas_beta_col), se_g = get(gwas_se_col),
             maf_g  = get(gwas_maf_col))],
    by = "snp_pos", all = FALSE
  )
  
  # 3) QC essentials
  m <- m[is.finite(beta_p) & is.finite(se_p) & is.finite(beta_g) & is.finite(se_g)]
  if (!nrow(m)) stop("No complete rows after QC.")
  
  # 4) Harmonize GWAS to pQTL ALT (effect allele)
  match_ok   <- m$p_ref == m$g_ref & m$p_alt == m$g_alt
  match_swap <- m$p_ref == m$g_alt & m$p_alt == m$g_ref
  
  if (any(match_swap)) {
    idx <- which(match_swap)
    m$beta_g[idx] <- -m$beta_g[idx]
    tmp <- m$g_ref[idx]; m$g_ref[idx] <- m$g_alt[idx]; m$g_alt[idx] <- tmp
    if (!all(is.na(m$maf_g[idx]))) m$maf_g[idx] <- 1 - m$maf_g[idx]
  }
  
  if (drop_ambiguous) {
    ambiguous <- is_ambiguous(m$p_ref, m$p_alt) | is_ambiguous(m$g_ref, m$g_alt)
  } else ambiguous <- rep(FALSE, nrow(m))
  
  keep <- (match_ok | match_swap) & !ambiguous
  m <- m[keep]
  if (!nrow(m)) stop("No aligned rows after harmonization.")
  
  # 5) Final unique key AFTER harmonization: chr:pos:ref:alt (pQTL alleles)
  m[, snpid := paste(chr, pos, p_ref, p_alt, sep=":")]
  
  # 6) Variances, MAF, p-values
  m[, varbeta_p := se_p^2]
  m[, varbeta_g := se_g^2]
  m[, MAF := fifelse(is.finite(maf_g), maf_g,
                     fifelse(is.finite(maf_p), maf_p, NA_real_))]
  m[, pval_p := 2*pnorm(-abs(beta_p/se_p))]
  m[, pval_g := 2*pnorm(-abs(beta_g/se_g))]
  
  # 7) Deduplicate exact variants (keep smallest variance)
  setorder(m, snpid, varbeta_p, varbeta_g)
  m <- m[, .SD[1], by = snpid]
  
  # 8) Build coloc datasets (ordered)
  setorder(m, chr, pos)
  d1 <- list( # pQTL
    snp      = m$snpid,
    beta     = m$beta_p,
    varbeta  = m$varbeta_p,
    MAF      = m$MAF,
    N        = N_pqtl,
    type     = "quant",
    position = m$pos
  )
  d2 <- list( # GWAS
    snp      = m$snpid,
    beta     = m$beta_g,
    varbeta  = m$varbeta_g,
    MAF      = m$MAF,
    N        = N_gwas,
    type     = "quant",
    position = m$pos
  )
  
  list(data1 = d1, data2 = d2, merged = m[])
}

extract_abf_summary <- function(abf) {
  s <- abf$summary
  if (is.null(s)) return(NULL)
  
  if (is.null(dim(s))) {
    # summary is a named vector
    list(
      nsnps = as.integer(s[["nsnps"]]),
      PP0   = as.numeric(s[["PP.H0.abf"]]),
      PP1   = as.numeric(s[["PP.H1.abf"]]),
      PP2   = as.numeric(s[["PP.H2.abf"]]),
      PP3   = as.numeric(s[["PP.H3.abf"]]),
      PP4   = as.numeric(s[["PP.H4.abf"]])
    )
  } else {
    # summary is a data frame with at least one row
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


###### RUN COLOC ######
# Define your protein loci and pQTL files
keep_pairs <- list(
  ASAT = c("THBS2"),
  VAT  = c("ABL1", "CCL17", "MSR1", "CFB", "KHK", "TSPAN8", "SELPLG",
           "ANXA2", "ASAH2", "ADAMTSL5", "SHBG", "ITGB6"),
  GFAT = c("SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC")
)

proteins <- fread("/Volumes/medpop_esp2/mpan/Projects/Adiposity/Adiposity_Omics/data/olink_protein_map_3k.tsv",
                  select = c("HGNC.symbol", "chr", "gene_start", "gene_end"))
proteins <- subset(proteins, proteins$HGNC.symbol %in% c("THBS2","ABL1", "CCL17", "MSR1", "CFB", "KHK", "TSPAN8", "SELPLG",
                            "ANXA2", "ASAH2", "ADAMTSL5", "SHBG", "ITGB6","SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC"))
colnames(proteins) <- c("gene", "chr", "start", "end")
proteins$pqtl_file <- paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/pqtl/",proteins$gene,"_cis_pqtl.txt.gz")

# Define your GWAS files
gwas_paths <- list(
  GFAT = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_gfatadjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz",
  ASAT = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_asatadjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz",
  VAT  = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_vatadjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz"
)

gwas_list <- lapply(gwas_paths, read_gwas)

# Sample sizes
gwas_N <- 32950

# Run the combination
results_dir <- "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/tissue_coloc/pqtl"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

out_rows <- list(); k <- 1L

for (i in seq_len(nrow(proteins))) {
  pr <- proteins[i]
  
  # get N directly from this protein's pQTL sumstats
  N_pqtl_i <- try(get_pqtl_N(pr$pqtl_file), silent = TRUE)
  if (inherits(N_pqtl_i, "try-error")) next
  
  for (depot in names(gwas_list)) {
    gwas_dt <- gwas_list[[depot]]
    
    inp <- try(prepare_coloc_inputs_pqtl_file(
      pqtl_path = pr$pqtl_file, gwas = gwas_dt,
      chr = pr$chr, region_start = pr$start, region_end = pr$end,
      N_pqtl = N_pqtl_i, N_gwas = gwas_N
    ), silent = TRUE)
    if (inherits(inp, "try-error")) next
    if (length(inp$data1$snp) < 5L) next
    
    abf <- try(coloc.abf(inp$data1, inp$data2), silent = TRUE)
    if (inherits(abf, "try-error")) next
    s <- extract_abf_summary(abf); if (is.null(s)) next
    
    m <- inp$merged
    out_rows[[k]] <- data.table(
      gene = pr$gene, depot = depot,
      chr = pr$chr, region_start = pr$start, region_end = pr$end,
      nsnps = s$nsnps,
      PP0 = s$PP0, PP1 = s$PP1, PP2 = s$PP2, PP3 = s$PP3, PP4 = s$PP4,
      min_p_pqtl = min(m$pval_p, na.rm = TRUE),
      min_p_gwas = min(m$pval_g, na.rm = TRUE),
      N_pqtl = N_pqtl_i,           
      N_gwas = gwas_N
    )
    k <- k + 1L
  }
}

res_pqtl <- data.table::rbindlist(out_rows, fill = TRUE)

# Filter res
res_filtered_pqtl <- res_pqtl[
  mapply(function(gwas_val, gene_val) {
    gwas_val %in% names(keep_pairs) &&
      gene_val %in% keep_pairs[[gwas_val]]
  }, depot, gene)
]

fwrite(res_pqtl, file.path(results_dir, "coloc_abf_pQTL-GWAS.2025-08-10.tsv"), sep = "\t")
fwrite(res_filtered_pqtl, file.path(results_dir, "coloc_abf_pQTL-GWAS_candidates.2025-08-10.tsv"), sep = "\t")


























###### VISUALIZE AND SUMMARIZE RESULTS ######
# res_filtered_pqtl <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/tissue_coloc/pqtl/coloc_abf_pQTL-GWAS_candidates.2025-08-10.tsv")
# 
# # Overall PP4 distribution
# summary(res_filtered_pqtl$PP4)
# hist(res_filtered_pqtl$PP4, breaks=50, main="PP4 distribution", xlab="PP4")
# 
# # Melt PP3 and PP4 into long format for plotting
# df_long <- melt(
#   res_filtered_pqtl,
#   id.vars = c("gene", "depot"),
#   measure.vars = c("PP3", "PP4"),
#   variable.name = "PP_type",
#   value.name = "PP_value"
# )
# 
# # Order genes alphabetically within each depot facet
# df_long <- df_long[order(depot, gene)]
# df_long$gene <- factor(df_long$gene, levels = unique(df_long$gene))
# 
# df_long <- df_long |>
#   dplyr::mutate(
#     PP_label = dplyr::if_else(PP_value < 0.01, "<0.01",
#                               sprintf("%.2f", round(PP_value, 2)))
#   )
# 
# svg("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure5b.svg",
#     width = 3.7, 
#     height = 5.5, 
#     pointsize = 7,    
#     family = "Arial")
# 
# # Ensure genes are ordered alphabetically within each depot facet
# # reorder_within() handles per-facet ordering; we don't need to mutate factor levels globally
# ggplot(
#   df_long,
#   aes(x = PP_value,
#       y = reorder_within(gene, gene, depot),
#       fill = PP_type)
# ) +
#   # Shading bands
#   annotate("rect", xmin = 0.6, xmax = 0.8,  ymin = -Inf, ymax = Inf,
#            fill = "khaki", alpha = 0.25, color = NA) +
#   annotate("rect", xmin = 0.8, xmax = 0.88, ymin = -Inf, ymax = Inf,
#            fill = "lightgreen", alpha = 0.25, color = NA) +
#   
#   # Bars
#   geom_col(position = position_dodge(width = 0.8), width = 0.7) +
#   
#   # Value labels placed just inside the panel edge if near 0.85
#   geom_text(
#     aes(x = pmin(PP_value + 0.01, 0.845), label = PP_label),
#     position = position_dodge(width = 0.8),
#     hjust = 0, size = 2.2, color = "black"
#   ) +
#   
#   facet_grid(rows = vars(depot), scales = "free_y", space = "free_y") +
#   scale_y_reordered(expand = expansion(add = c(0.6, 0.6), mult = c(0, 0))) +
#   
#   scale_x_continuous(
#     limits = c(0, 0.88),
#     expand = expansion(mult = 0, add = c(0.01, 0))
#   ) +
#   
#   scale_fill_manual(values = c("PP3" = "lightgrey", "PP4" = "#4AA1A1")) +
#   labs(x = "Posterior Probability", y = "", fill = "Hypothesis",
#        title = "PP3 vs. PP4 from pQTL colocalization analysis") +
#   theme_bw(base_size = 7) +
#   theme(
#     axis.text.x  = element_text(size = 6),
#     axis.text.y  = element_text(size = 6),
#     axis.title.x = element_text(size = 7),
#     axis.title.y = element_blank(),
#     strip.text   = element_text(size = 7, margin = margin(1, 1, 1, 1)),
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
# 
# 
# dev.off()
# 
# # a smaller version of the graph with only proteins of interest (those with moderate to strong evidence of either PP3 or PP4)
# small <- subset(df_long, df_long$gene %in% c("SHBG", "CFB"))
# small <- subset(small, !(small$gene=="SHBG" & small$depot=="VAT"))
# 
# svg("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/Figure5b_small.svg",
#     width = 3.5, 
#     height = 2, 
#     pointsize = 7,    
#     family = "Arial")
# 
# ggplot(
#   small,
#   aes(x = PP_value,
#       y = reorder_within(gene, gene, depot),
#       fill = PP_type)
# ) +
#   # Shading bands
#   annotate("rect", xmin = 0.6, xmax = 0.8,  ymin = -Inf, ymax = Inf,
#            fill = "khaki", alpha = 0.25, color = NA) +
#   annotate("rect", xmin = 0.8, xmax = 0.88, ymin = -Inf, ymax = Inf,
#            fill = "lightgreen", alpha = 0.25, color = NA) +
#   
#   # Bars
#   geom_col(position = position_dodge(width = 0.8), width = 0.7) +
#   
#   # Value labels placed just inside the panel edge if near 0.85
#   geom_text(
#     aes(x = pmin(PP_value + 0.01, 0.845), label = PP_label),
#     position = position_dodge(width = 0.8),
#     hjust = 0, size = 2.2, color = "black"
#   ) +
#   
#   facet_grid(rows = vars(depot), scales = "free_y", space = "free_y") +
#   scale_y_reordered(expand = expansion(add = c(0.6, 0.6), mult = c(0, 0))) +
#   
#   scale_x_continuous(
#     limits = c(0, 0.88),
#     expand = expansion(mult = 0, add = c(0.01, 0))
#   ) +
#   
#   scale_fill_manual(values = c("PP3" = "lightgrey", "PP4" = "#4AA1A1")) +
#   labs(x = "Posterior Probability", y = "", fill = "Hypothesis",
#        title = "PP3 vs. PP4 from pQTL colocalization analysis") +
#   theme_bw(base_size = 7) +
#   theme(
#     axis.text.x  = element_text(size = 6),
#     axis.text.y  = element_text(size = 6),
#     axis.title.x = element_text(size = 7),
#     axis.title.y = element_blank(),
#     strip.text   = element_text(size = 7, margin = margin(1, 1, 1, 1)),
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
#   
# dev.off()
