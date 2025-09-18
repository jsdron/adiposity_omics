#########################################
# Script: locuszoom_NFASC.R
# Description: Generates multi-panel locuszoom plots for the NFASC locus,
#              integrating GWAS, eQTL, and pQTL data with LD information 
#              and gene tracks.
#
# Key Outputs:
#   - Figure (PDF; locuszoom plots of NFASC across GWAS, eQTL, pQTL, and gene tracks)
#########################################

###### LIBRARIES ######
library(data.table)
library(dplyr)
library(locuszoomr)
library(AnnotationHub)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ensembldb")
#BiocManager::install("EnsDb.Hsapiens.v86")

library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86

###### VARIABLES ######
gene <- "NFASC"
tissue_type <- "Brain_Nucleus_accumbens_basal_ganglia"
fat_depot <- "gfat"
#index_snp <- "rs2955617" 
ld_token <- "fbe14f7cef7a"

###### PREPARE DATA ######
# eQTL -----
eqtl_main <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/target_protein.all_tissues_eqtls.tsv.gz")

eqtl_df <- subset(eqtl_main, eqtl_main$gene_symbol==gene & eqtl_main$tissue==tissue_type)
eqtl_df <- eqtl_df[,c("chr", "pos", "ref", "alt", "af", "slope", "slope_se", "pval_nominal", "ma_samples")]
colnames(eqtl_df) <- c("chr", "pos", "NEA", "EA", "EAF_eQTL", "b_eQTL", "se_eQTL", "pval_eQTL", "N_eQTL")
eqtl_df$chr <- gsub("chr", "", eqtl_df$chr)
eqtl_df$chr <- as.integer(eqtl_df$chr)

setDT(eqtl_df)

eqtl_df[, `:=`(NEA = toupper(NEA), EA = toupper(EA))]

# Alphabetize alleles and realign beta/EAF to the (new) EA
eqtl_df[, `:=`(A1 = pmin(NEA, EA), A2 = pmax(NEA, EA))]
eqtl_df[, swapped := EA != A2 & !is.na(EA) & !is.na(A2)]
eqtl_df[swapped == TRUE, `:=`(b_eQTL = -b_eQTL, EAF_eQTL = 1 - EAF_eQTL)]

eqtl_df[, `:=`(NEA = A1, EA = A2)]
eqtl_df[, c("A1","A2","swapped") := NULL]

# GWAS -----
gwas_df <- fread(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/GWAS/bolt_lmm/for_pub/GWAS_",fat_depot,"adjbmi3_UKB-noProt_EUR_hg38_INFOgt0.3_MAFgt0.005.stats.gz"),
                 select = c("SNP", "CHR", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "P_BOLT_LMM"))
colnames(gwas_df) <- c("SNP","chr", "pos", "NEA", "EA", "EAF_outcome", "b_outcome", "se_outcome", "pval_outcome")
gwas_df$N_outcome <- 32950 # the same size for all the fat depot GWAS

setDT(gwas_df)

gwas_df[, `:=`(NEA = toupper(NEA), EA = toupper(EA))]

# Alphabetize alleles and realign beta/EAF to the (new) EA
gwas_df[, `:=`(A1 = pmin(NEA, EA), A2 = pmax(NEA, EA))]
gwas_df[, swapped := EA != A2 & !is.na(EA) & !is.na(A2)]
gwas_df[swapped == TRUE, `:=`(b_outcome = -b_outcome, EAF_outcome = 1 - EAF_outcome)]

gwas_df[, `:=`(NEA = A1, EA = A2)]
gwas_df[, c("A1","A2","swapped") := NULL]

# pQTL -----
pqtl_df <- fread(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/pqtl/",gene,"_cis_pqtl.txt.gz"),
                 select = c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P", "N"))
colnames(pqtl_df) <- c("chr", "pos", "NEA", "EA", "EAF_pQTL", "b_pQTL", "se_pQTL", "pval_pQTL", "N_pQTL")
pqtl_df$pval_pQTL <- 10^(-pqtl_df$pval_pQTL)

setDT(pqtl_df)

pqtl_df[, `:=`(NEA = toupper(NEA), EA = toupper(EA))]

# Alphabetize alleles and realign beta/EAF to the (new) EA
pqtl_df[, `:=`(A1 = pmin(NEA, EA), A2 = pmax(NEA, EA))]
pqtl_df[, swapped := EA != A2 & !is.na(EA) & !is.na(A2)]
pqtl_df[swapped == TRUE, `:=`(b_pQTL = -b_pQTL, EAF_pQTL = 1 - EAF_pQTL)]

pqtl_df[, `:=`(NEA = A1, EA = A2)]
pqtl_df[, c("A1","A2","swapped") := NULL]

###### MERGE DATA ######
merged_df <- merge(gwas_df, eqtl_df, by=c("chr", "pos", "NEA", "EA"), all = FALSE)
merged_df <- merge(merged_df, pqtl_df, by=c("chr", "pos", "NEA", "EA"), all = FALSE)

###### PLOT PREP ######
# eQTL -----
dat.eqtl <- data.frame(SNP = merged_df$SNP, 
                       eaf = merged_df$EAF_eQTL,
                       chrom = merged_df$chr,
                       pos = merged_df$pos,
                       LOG10P = -log10(merged_df$pval_eQTL), 
                       pval = merged_df$pval_eQTL
)
dat.eqtl <- subset(dat.eqtl, is.finite(pval) & pval > 0)

loc.eqtl <- locus(data = dat.eqtl, 
                  gene = gene, 
                  chrom = "chrom", 
                  pos = "pos", 
                  xrange = NULL, seqname = NULL, flank = NULL,
                  ens_db = ensdb, 
                  fix_window = 1e6, 
                  labs ="SNP", 
                  p = "pval", 
                  #index_snp = index_snp,
                  std_filter = TRUE
)

loc.eqtl <- link_LD(loc.eqtl, 
                    pop = c("CEU", "GBR", "TSI", "FIN", "IBS"),
                    genome_build = "GRCh38",
                    token = ld_token)
# loc.eqtl$data  <- subset(loc.eqtl$data,  SNP == loc.eqtl$index | !is.na(ld)) # Remove NAs

# GWAS -----
dat.gwas <- data.frame(SNP = merged_df$SNP, 
                       eaf = merged_df$EAF_outcome,
                       chrom = merged_df$chr,
                       pos = merged_df$pos,
                       LOG10P = -log10(merged_df$pval_outcome), 
                       pval = merged_df$pval_outcome
)
dat.gwas <- subset(dat.gwas, is.finite(pval) & pval > 0)

loc.gwas <- locus(data = dat.gwas, 
                  gene = gene, 
                  chrom = "chrom", 
                  pos = "pos", 
                  xrange = NULL, seqname = NULL, flank = NULL,
                  ens_db = ensdb, 
                  fix_window = 1e6, 
                  labs ="SNP", 
                  p = "pval", 
                  #index_snp = index_snp,
                  std_filter = TRUE
)

loc.gwas <- link_LD(loc.gwas, 
                    pop = c("CEU", "GBR", "TSI", "FIN", "IBS"),
                    genome_build = "GRCh38",
                    token = ld_token)

# pQTL -----
dat.pqtl <- data.frame(SNP = merged_df$SNP, 
                       eaf = merged_df$EAF_pQTL,
                       chrom = merged_df$chr,
                       pos = merged_df$pos,
                       LOG10P = log10(merged_df$pval_pQTL), 
                       pval = merged_df$pval_pQTL
)
dat.pqtl <- subset(dat.pqtl, is.finite(pval) & pval > 0)

loc.pqtl <- locus(data = dat.pqtl, 
                  gene = gene, 
                  chrom = "chrom", 
                  pos = "pos", 
                  xrange = NULL, seqname = NULL, flank = NULL,
                  ens_db = ensdb, 
                  fix_window = 1e6, 
                  labs ="SNP", 
                  p = "pval", 
                  #index_snp = index_snp,
                  std_filter = TRUE
)

loc.pqtl <- link_LD(loc.pqtl, 
                    pop = c("CEU", "GBR", "TSI", "FIN", "IBS"),
                    genome_build = "GRCh38",
                    token = ld_token)


###### PLOT ######
ld_cols <- c("grey80","navyblue","turquoise3",
             "forestgreen","orange","red3","purple3")
ld_labels <- c("NA","0.0 - 0.2","0.2 - 0.4",
               "0.4 - 0.6","0.6 - 0.8","0.8 - 1.0","Index SNP")

p1 <- gg_scatter(loc.gwas, size = 1, index_snp = "index_snp", labels = "index", xticks = FALSE, legend_pos = "topleft") + 
  labs(y = expression(GWAS ~ -log[10](P))) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank(),
        legend.margin        = margin(0,0,0,0),
        legend.box.margin    = margin(0,0,0,0),
        legend.key = element_rect(fill="transparent", colour=NA),
        legend.box.background = element_rect(fill="transparent", colour="black"),
        legend.key.height = grid::unit(6, "pt"),
        legend.key.width  = grid::unit(7, "pt"),
        legend.spacing.y  = grid::unit(0, "pt"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.25, colour = "black"), # ticks match axis 
        axis.line = element_line(linewidth = 0.25, colour = "black"), # thinner axis lines
        legend.text = element_text(size = 6)
  ) + 
  scale_fill_manual(values=ld_cols, labels=ld_labels, drop=FALSE, name=expression(LD~r^2))

p2 <- gg_scatter(loc.eqtl, size = 1, index_snp = "index_snp", labels = NULL, xticks = FALSE, legend_pos = NULL) + 
  labs(y = expression(eQTL ~ -log[10](P))) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values=ld_cols, labels=ld_labels, drop=FALSE, name=expression(LD~r^2)) +
  theme(axis.title = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.25, colour = "black"), # ticks match axis 
        axis.line = element_line(linewidth = 0.25, colour = "black"), # thinner axis lines
        axis.text = element_text(size = 6)
  )

p3 <- gg_scatter(loc.pqtl, size = 1, index_snp = "index_snp", labels = NULL, legend_pos = NULL) + 
  labs(y = expression(pQTL ~ -log[10](P))) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values=ld_cols, labels=ld_labels, drop=FALSE, name=expression(LD~r^2)) +
  theme(axis.title = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.25, colour = "black"), # ticks match axis 
        axis.line = element_line(linewidth = 0.25, colour = "black"), # thinner axis lines
        axis.text = element_text(size = 6)
  )


p5 <- gg_genetracks(cex.axis = 0, cex.lab = 0, cex.text = 0.4,
                    loc.eqtl, xticks = FALSE,
                    filter_gene_biotype = "protein_coding", maxrows = 15, 
                    highlight = gene,
                    italics = TRUE, gene_col = "grey", exon_col = "darkgrey", exon_border = "black"
)

ap_full <- cowplot::align_plots(
  p1 + theme(plot.margin = margin(5.5,5.5,5.5,5.5)),
  p2 + theme(plot.margin = margin(5.5,5.5,5.5,5.5)),
  p3 + theme(plot.margin = margin(5.5,5.5,5.5,5.5)),
  p5 + theme(plot.margin = margin(5.5,5.5,5.5,5.5)),
  align = "v", axis = "l"
)

# cowplot::plot_grid(ap_full[[1]], ap_full[[2]], ap_full[[3]], ncol = 1, rel_heights = c(1,1,1.05))
cowplot::plot_grid(ap_full[[1]], ap_full[[2]], ap_full[[3]], ap_full[[4]], ncol = 1,  rel_heights = c(1,1,1.05,0.6))

###### SAVE FIGURE ######
pdf(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/manuscript/pdf/Figure5_", gene, "-locuszoom.pdf"), 
    width = 3.5, 
    height = 7, 
    family = "Arial") 
cowplot::plot_grid(ap_full[[1]], ap_full[[2]], ap_full[[3]], ap_full[[4]], ncol = 1,  rel_heights = c(1,1,1.05,1.8))
dev.off()

