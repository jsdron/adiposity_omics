#########################################
# Script: Out-Fat_Exp-Metab.R
# Description: Two-sample Mendelian randomization (MR) pipeline for metabolite exposures.
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

setDTthreads(0)

###### SET VARIABLES ######
path <- "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics"
run <- "MET_ADIPOSE_twosample/"
date <- "2025-02-08_test"

###### LOG FILE SETUP ######
log_file <- paste0(path, "/src/log/MET_Adipose_twosample_", date)
log_connection <- file(log_file, open = "a")
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  writeLines(paste(timestamp, message), log_connection)
}

###### LOAD OUTCOME AND EXPOSURE DEFINITIONS ######
outcome_file <- fread(file = paste0(path, "/data/Met_Adipose_outcomes_twosample.csv"))
outcome_list <- outcome_file$custom_name

exposure_file <- fread(file = paste0(path, "/data/raw_met-d_mp_samplesize.txt"))
exposure_file <- subset(exposure_file, unit != "ratio")
exposure_list <- exposure_file$metabolite
exposure_names <- exposure_file$clean_name

all_outcome_results <- NULL

###### MR ANALYSIS LOOP ######
# Run MR, sensitivity analysis, and horizontal pleiotropy test for CAD outcome with looping structure intact
for (j in 1:length(outcome_list)) {
  print("Loading summary statistics...")
  log_message("Loading summary statistics...")
  outcome_sumstats <- fread(outcome_file$path[j], header = TRUE, stringsAsFactors = FALSE)
  print("Formatting outcome data...")
  log_message("Formatting outcome data...")
  outcome_sumstats$phenotype <- outcome_list[j]
  outcome_sumstats <- data.frame(outcome_sumstats)
  outcome <- format_data(outcome_sumstats, type = "outcome",
                        snp_col = "SNP",
                        beta_col = "BETA",
                        se_col = "SE",
                        eaf_col = "A1FREQ",
                        effect_allele_col = "ALLELE1",
                        other_allele_col = "ALLELE0",
                        pval_col = "P_BOLT_LMM",
                        phenotype_col = "phenotype")
  
  mr_output <- NULL
  sensitivity_output <- NULL
  pleiotropy_output <- NULL
  all_dat <- NULL
  
  # Creating dataframe for instrument lists
  instruments_final = data.frame(1, "rs001")
  colnames(instruments_final) = c("id", "filler")
  
  for (i in 1:length(exposure_list)) {
    print(paste0(i, ": ", exposure_names[i], " for ", j, ": ", outcome_list[j]))
    log_message(paste0(i, ": ", exposure_names[i], " for ", j, ": ", outcome_list[j]))
    
    ### MR: Exposure -> Outcome ###
    exposure <- extract_instruments(exposure_list[i]) # they are GWS loci
    dat <- harmonise_data(exposure, outcome)
    colnames(dat)[colnames(dat) %in% c("id.outcome")] <- c("outcome_alt")
    dat$id.outcome <- dat$outcome
    
    dat$samplesize.exposure = exposure_file$samplesize[i]
    dat$samplesize.outcome = outcome_file$samplesize[j]
    dat = steiger_filtering(dat)
    dat = dat[which(dat$steiger_dir == TRUE & dat$steiger_pval < 0.05),]
    
    dat <- dat[order(dat$pval.exposure),] # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
    dat <- dat[!duplicated(dat$SNP),] # Ordering the table by pval first ensures that the more significant SNP is kept
    dat <- dat[!duplicated(dat$pos.exposure),]
    dat <- clump_data(dat, clump_p1 = 5e-8, clump_p2 = 5e-8) # defaults here: https://mrcieu.github.io/TwoSampleMR/reference/clump_data.html

    all_dat <- rbind(all_dat, dat) # all the data pre-MR
    
    res <- mr(dat) 
    res <- generate_odds_ratios(res)
    res$exposure <- res$id.exposure
  
    # MR Raps method
    res2 <- mr.raps(dat$beta.exposure[dat$mr_keep==T], dat$beta.outcome[dat$mr_keep==T], dat$se.exposure[dat$mr_keep==T], dat$se.outcome[dat$mr_keep==T])
    
    res2 <- data.frame(res2)
    res2$id.exposure <- res$id.exposure[1]
    res2$id.outcome <- outcome_list[j]
    res2$exposure <- res$id.exposure[1]
    res2$outcome <- outcome_list[j]
    res2$method <- "MR Raps"
    res2$nsnp <- length(dat$beta.exposure[dat$mr_keep==T])
    colnames(res2)[colnames(res2) %in% c("beta.hat", "beta.se", "beta.p.value")] <- c("b", "se", "pval")
    res2 <- res2[c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval")]
    res2 <- generate_odds_ratios(res2)
    
    # Aggregating results
    mr_output <- rbind(mr_output, res)
    mr_output <- rbind(mr_output, res2) 
    
    ### Sensitivity analysis ###
    sens <- mr_heterogeneity(dat)
    sensitivity_output <- rbind(sensitivity_output, sens) # all the sensitivity results
    
    ### Horizontal pleiotropy test ###
    pleio <- mr_pleiotropy_test(dat)
    pleiotropy_output <- rbind(pleiotropy_output, pleio) # all the pleiotropy results
    
    p1 <- mr_scatter_plot(res, dat)
    ggsave(p1[[1]], file <- paste0(path, "/results/MR/", run, outcome_list[j],"/scatter_plot/",outcome_list[j],"_",exposure_names[i],".pdf"), width = 7, height = 7)
    
    #Appending the instruments for each run now
    snp_append <- unique(subset(dat, mr_keep == TRUE)$SNP)
    exposure_new <- data.frame(1:length(snp_append), snp_append)
    exposure_clean_name <- gsub("met-d-", "", exposure_list[i])
    colnames(exposure_new) <- c("id", exposure_clean_name)
    instruments_final <- merge(instruments_final, exposure_new, by = "id", all = TRUE)

  }

  mr_output <- as.data.frame(mr_output)
  sensitivity_output <- as.data.frame(sensitivity_output)
  pleiotropy_output <- as.data.frame(pleiotropy_output)
  
  # Writing out files: 1 per outcome
  fwrite(mr_output, paste0(path, "/results/MR/", run, outcome_list[j], "/mr_results_", outcome_list[j], "_", date, ".csv"))
  fwrite(sensitivity_output, paste0(path, "/results/MR/", run, outcome_list[j], "/mr_sensitivity_", outcome_list[j], "_", date, ".csv"))
  fwrite(pleiotropy_output, paste0(path, "/results/MR/", run, outcome_list[j], "/mr_pleiotropy_", outcome_list[j], "_", date, ".csv"))
  
  instruments_final <- as.data.frame(instruments_final)
  instruments_final <- instruments_final[,-2]
  fwrite(instruments_final, paste0(path, "/results/MR/", run, "genetic_instruments_",outcome_list[j], "_", date, ".csv")) 
  
  all_outcome_results <- rbind(all_outcome_results, mr_output)
  
  outcome_sumstats <- NULL
  outcome <- NULL
  
}

fwrite(all_outcome_results, paste0(path, "/results/MR/", run, "all_outcomes_mr_results_Met-Adipose", date, ".csv")) 

close(log_connection)

print("DONE")
