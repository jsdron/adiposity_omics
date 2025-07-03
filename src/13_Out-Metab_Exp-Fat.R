#########################################
# Script: Out-Metab_Exp-Fat.R
# Description: Two-sample Mendelian randomization (MR) pipeline for fat exposures on metabolites.
# Key Outputs:
#   - MR results per outcome (CSV)
#   - Sensitivity analyses (heterogeneity, pleiotropy)
#   - Genetic instruments used
#   - Summary file of all outcome MR results
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(stringr)
library(ggplot2)
library(R.utils)
library(ieugwasr)
library(MRInstruments)
library(mr.raps)
library(vcfR)

setDTthreads(0)

###### SET VARIABLES ######
path <- "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics"
run <- "ADIPOSE_MET_twosample/"
date <- "2025-03-16"

###### LOG FILE SETUP ######
log_file <- paste0(path, "/src/log/Adipose_MET_twosample_", date)
log_connection <- file(log_file, open = "a")
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  writeLines(paste(timestamp, message), log_connection)
}

###### LOAD OUTCOME AND EXPOSURE DEFINITIONS ######
# This file only has the exposure of interest.
exposure_file <- fread(file = paste0(path, "/data/Met_Adipose_outcomes_twosample.csv"), header = TRUE, stringsAsFactors = FALSE)
exposure_list <- exposure_file$custom_name   

# Listed based on forward MR
met_list <- c("M_HDL_CE", "HDL_P", "M_HDL_FC", "M_HDL_C", "ApoA1", "M_HDL_PL", "M_HDL_P", "M_HDL_L", "HDL_PL", "HDL_L", "HDL_CE")

outcome_file <- fread(file = paste0(path, "/data/raw_met-d_mp_samplesize.txt"), header = TRUE, stringsAsFactors = FALSE)
outcome_file <- subset(outcome_file, clean_name %in% met_list)

chr_list <- list('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')

# Sub-setting out percentage values
outcome_file <- subset(outcome_file, outcome_file$unit!="ratio")
outcome_list = outcome_file$metabolite
outcome_names = outcome_file$clean_name

outcome_path_list <- NULL
for (i in 1:length(outcome_list)){
  outcome_path_list[i] <- paste0("/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/data/", outcome_list[i], ".vcf.gz")
}

all_outcome_results <- NULL

###### MR ANALYSIS LOOP ######
i <- 1
j <- 1

# Run MR, sensitivity analysis, and horizontal pleiotropy test for CAD outcome with looping structure intact
for (j in 1:length(exposure_list)) {
  print("Loading summary statistics...")
  log_message("Loading summary statistics...")
  exposure_sumstats <- fread(exposure_file$path[j], header = TRUE, stringsAsFactors = FALSE)
  print("Formatting exposure data...")
  log_message("Formatting exposure data...")
  exposure_sumstats$phenotype <- exposure_list[j]
  exposure_sumstats$CHR <- ifelse(exposure_sumstats$CHR == 23, "X", exposure_sumstats$CHR)
  
  # Processing adipose sum stats to filter by P value for significance - ONLY 34 LEFT
  exposure_sumstats <- exposure_sumstats[exposure_sumstats$P_BOLT_LMM < 5e-8]

  exposure_sumstats <- data.frame(exposure_sumstats)
  exposure <- format_data(exposure_sumstats, type = "exposure",
                        snp_col = "SNP",
                        beta_col = "BETA",
                        se_col = "SE",
                        eaf_col = "A1FREQ",
                        effect_allele_col = "ALLELE1",
                        other_allele_col = "ALLELE0",
                        pval_col = "P_BOLT_LMM",
                        phenotype_col = "phenotype",
                        chr_col="CHR",
                        pos_col="BP")
  
  colnames(exposure)[colnames(exposure) %in% c("id.exposure")] <- c("exposure_alt")
  exposure$id.exposure <- exposure$exposure
  
  mr_output <- NULL
  sensitivity_output <- NULL
  pleiotropy_output <- NULL
  all_dat <- NULL
  
  # Creating dataframe for instrument lists
  instruments_final <- data.frame(1, "rs001")
  colnames(instruments_final) <- c("id", "filler")
  
  for (i in 1:length(outcome_list)) {
    print(paste0(i, ": ", outcome_names[i], " for ", j, ": ", exposure_list[j]))
    log_message(paste0(i, ": ", outcome_names[i], " for ", j, ": ", exposure_list[j]))
    
    # Downloading metabolite sum statss
    vcf <- read.vcfR(outcome_path_list[i], verbose = FALSE)
    vcf_df <- as.data.frame(vcf@fix)
    head(vcf_df)
    gt_df <- as.data.frame(vcf@gt)
    head(gt_df)
    vcf_full_df <- cbind(vcf_df, gt_df)
    head(vcf_full_df)
    
    # Filtering to the overlapping positions just to save time
    vcf_overlap <- vcf_full_df[vcf_full_df$POS %in% exposure$pos.exposure,]  
    
    # Split column into five separate columns
    vcf_overlap <- vcf_overlap %>%
      separate(outcome_list[i], into = c("ES", "SE", "LP", "AF", "ID"), sep = ":")
    
    vcf_overlap$CHROM <- ifelse(vcf_overlap$CHROM == 23, "X", vcf_overlap$CHROM) # Transferring 23rd chromosome to "X" for consistency
    vcf_overlap <- unique(vcf_overlap)
    
    if (is.null(vcf_overlap) || nrow(vcf_overlap) == 0) {
      print(paste0("Skipping ", outcome_names[i], ": vcf_overlap null"))
      log_message(paste0("Skipping ", outcome_names[i], ": vcf_overlap null"))
    } else{
      exposure <- exposure[exposure$pos.exposure %in% vcf_overlap$POS,] # Only selecting the variants that are overlapping
      exposure <- unique(exposure)
      
      if (is.null(exposure) || nrow(exposure) == 0) {
        print(paste0("Skipping ", outcome_names[i], ": exposure overlap null"))
        log_message(paste0("Skipping ", outcome_names[i], ": exposure overlap null"))
      } else{
        vcf_overlap$phen <- outcome_names[i]
        vcf_overlap <- data.frame(vcf_overlap)
        vcf_overlap$LP <- as.numeric(vcf_overlap$LP)
        outcome <- format_data(vcf_overlap, type="outcome", 
                                     phenotype_col="phen", 
                                     snp_col="ID", 
                                     beta_col="ES", 
                                     se_col="SE", 
                                     eaf_col="AF",
                                     effect_allele_col="ALT", 
                                     other_allele_col="REF", 
                                     pval_col="LP", 
                                     chr_col="CHROM", 
                                     pos_col="POS",
                                     log_pval = T)
        colnames(outcome)[colnames(outcome) %in% c("id.outcome")] <- c("outcome_alt")
        outcome$id.outcome <- outcome$outcome
        outcome$chr.outcome <- as.integer(outcome$chr.outcome)
        outcome$pos.outcome <- as.integer(outcome$pos.outcome)
        rm(vcf, vcf_df, vcf_full_df, gt_df, vcf_overlap)
    
    dat <- harmonise_data(exposure, outcome)
    
    if (is.null(dat) || nrow(dat) == 0){
      print(paste0("Skipping ", outcome_list[i], ": Nothing after harmonize"))
      log_message(paste0("Skipping ", outcome_list[i], ": Nothing after harmonize"))
      next
    }
    
    dat$samplesize.exposure <- exposure_file$samplesize[j]
    dat$samplesize.outcome <- outcome_file$samplesize[i]
    dat <- steiger_filtering(dat)
    dat <- dat[which(dat$steiger_dir == TRUE & dat$steiger_pval < 0.05),]
    
    if (is.null(dat) || nrow(dat) == 0){
      print(paste0("Skipping ", outcome_list[i], ": Nothing passed steiger filtering"))
      log_message(paste0("Skipping ", outcome_list[i], ": Nothing passed steiger filtering"))
      next
    }
    
    dat <- dat[order(dat$pval.exposure),] # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
    dat <- dat[!duplicated(dat$SNP),] # Ordering the table by pval first ensures that the more significant SNP is kept
 
    clump <- ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure), # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = 0.001, clump_p = 5e-8,
                      bfile = "/Users/michael.tian/Desktop/Natarajan_Lab/Tools/g1000_eur")
    dat <- dat[dat$SNP %in% clump$rsid,]
    
    all_dat <- rbind(all_dat, dat) # all the data pre-MR
    
    res <- mr(dat) 
    res <- generate_odds_ratios(res)
    res$exposure <- res$id.exposure
    
    # MR Raps method
    res2 <- mr.raps(dat$beta.exposure[dat$mr_keep == T], dat$beta.outcome[dat$mr_keep == T], dat$se.exposure[dat$mr_keep == T], dat$se.outcome[dat$mr_keep == T])
    
    res2 <- data.frame(res2)
    res2$id.exposure <- res$id.exposure[1]
    res2$id.outcome <- outcome_names[j]
    res2$exposure <- res$id.exposure[1]
    res2$outcome <- outcome_names[j]
    res2$method <- "MR Raps"
    res2$nsnp <- length(dat$beta.exposure[dat$mr_keep == T])
    colnames(res2)[colnames(res2) %in% c("beta.hat", "beta.se", "beta.p.value")] <- c("b", "se", "pval")
    res2 <- res2[c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval")]
    res2 <- generate_odds_ratios(res2)
    
    # Aggregating results
    mr_output <- rbind(mr_output, res)
    mr_output <- rbind(mr_output, res2) 
    
    # Sensitivity analysis 
    sens <- mr_heterogeneity(dat)
    sensitivity_output <- rbind(sensitivity_output, sens) # all the sensitivity results
    
    # Horizontal pleiotropy test 
    pleio <- mr_pleiotropy_test(dat)
    pleiotropy_output <- rbind(pleiotropy_output, pleio) # all the pleiotropy results
    
    p1 <- mr_scatter_plot(res, dat)
    ggsave(p1[[1]], file <- paste0(path, "/results/MR/", run, exposure_list[j], "/scatter_plot/", exposure_list[j], "_", outcome_names[i], ".pdf"), width = 7, height = 7)
    
    #Appending the instruments for each run now
    snp_append <- unique(subset(dat, mr_keep == TRUE)$SNP)
    exposure_new <- data.frame(1:length(snp_append), snp_append)
    exposure_clean_name <- gsub("met-d-", "", outcome_list[i])
    colnames(exposure_new) <- c("id", exposure_clean_name)
    instruments_final <- merge(instruments_final, exposure_new, by = "id", all = TRUE)
    
  }
  mr_output <- as.data.frame(mr_output)
  sensitivity_output <- as.data.frame(sensitivity_output)
  pleiotropy_output <- as.data.frame(pleiotropy_output)
  
  # Writing out files: 1 per outcome
  fwrite(mr_output, paste0(path, "/results/MR/", run, exposure_list[j], "/mr_results_", exposure_list[j], "_", date, ".csv")) 
  fwrite(sensitivity_output, paste0(path, "/results/MR/", run, exposure_list[j], "/mr_sensitivity_", exposure_list[j], "_", date, ".csv")) 
  fwrite(pleiotropy_output, paste0(path, "/results/MR/", run, exposure_list[j], "/mr_pleiotropy_", exposure_list[j], "_", date, ".csv")) 
  
  instruments_final <- as.data.frame(instruments_final)
  instruments_final <- instruments_final[,-2]
  fwrite(instruments_final, paste0(path, "/results/MR/", run, "genetic_instruments_", exposure_list[j], "_", date, ".csv")) 
  
  all_outcome_results <- rbind(all_outcome_results, mr_output)
  
  outcome_sumstats <- NULL
  outcome <- NULL
  
    }
  }
}
test <- unique(all_outcome_results)

fwrite(all_outcome_results, paste0(path, "/results/MR/", run, "all_outcomes_mr_results_Adipose-MET", date, ".csv")) 

close(log_connection)

print("DONE")
