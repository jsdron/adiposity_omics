#########################################
# Script: Out-Prot_Exp-Fat.R
# Description: Two-sample Mendelian randomization (MR) pipeline for fat exposures on proteins
# Key Outputs:
#   - MR results per outcome (CSV)
#   - Sensitivity analyses (heterogeneity, pleiotropy)
#   - Genetic instruments used
#   - Summary file of all outcome MR results
#########################################

###### LIBRARIES ######
library(data.table)
library(synapser)
library(synapserutils)
library(TwoSampleMR)
library(tidyverse)
library(ggplot2)
library(R.utils)
library(MRInstruments)
library(tidyr)
library(dplyr)
library(ieugwasr)
library(genetics.binaRies)
library(MendelianRandomization)
library(mr.raps)

###### SET VARIABLES ######
path <- "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics"
run <- "ADIPOSE_PROT_twosample/"
date <- "2025-04-23"

###### LOG FILE SETUP ######
log_file <- paste0(path, "/src/log/Adipose_Prot_TwoSample_", date)
log_connection <- file(log_file, open = "a")
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  writeLines(paste(timestamp, message), log_connection)
}

###### SYNAPSE ######
synLogin(authToken = "XXX") # Replace 'XXX' with User's authToken fron Synapse
ref_rsid <- fread(paste0(path, "/data/hg38_common_chrpos_X.txt")) # This is a linker file to match RSIDs with chromosome positions for the X chromosome
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

###### LOAD OUTCOME AND EXPOSURE DEFINITIONS ######
sumstats_info <- fread(paste0(path, "/data/olink_protein_map_3k.tsv")) # Information on the sum stats / protein codes etc.; functions as a linker file
sumstats_info[1014, 5] <- "DKK1"
sumstats_info[1014, 9] <- "DKK1"
sumstats_info[414, 5] <- "SIT1"
sumstats_info[414, 9] <- "SIT1"
sumstats_info[215, 5] <- "DKK3"
sumstats_info[215, 9] <- "DKK3"
sumstats_info[882, 5] <- "DKK4"
sumstats_info[882, 9] <- "DKK4"
sumstats_info[1157, 5] <- "LSM1"
sumstats_info[1157, 9] <- "LSM1"
sumstats_info[1255, 5] <- "VNN2"
sumstats_info[1255, 9] <- "VNN2"
sumstats_info[1547, 5] <- "NPR1"
sumstats_info[1547, 9] <- "NPR1"
sumstats_info[2142, 5] <- "VNN1"
sumstats_info[2142, 9] <- "VNN1"
sumstats_info[2701, 5] <- "LSM8"
sumstats_info[2701, 9] <- "LSM8"
sumstats_info[2824, 5] <- "ERN1"
sumstats_info[2824, 9] <- "ERN1"
sumstats_info[2881, 5] <- "TOP1"
sumstats_info[2881, 9] <- "TOP1"

targets <- c("LEP", "ADIPOQ") # Set the proteins to run through MR batch here
sumstats_info <- sumstats_info[sumstats_info$Assay %in% targets, ]
outcome_names <- sumstats_info$Assay

exposure_file <- fread(file = paste0(path, "/data/Prot_Adipose_outcomes_twosample.csv"), header = TRUE, stringsAsFactors = FALSE)
exposure_list <- exposure_file$custom_name  

all_outcome_results <- NULL

###### DOWNLOAD SUMMARY STATS ######
l <- 1
syn_code_path_list <- list()
syn_code_cache_list <- list()

for (k in sumstats_info$Code){
  print(paste0("Downloading exposure: ", sumstats_info[sumstats_info$Code == k, ]$Assay[1]))
  log_message(paste0("Downloading exposure: ", sumstats_info[sumstats_info$Code == k, ]$Assay[1]))
  syn_code <- synGet(entity=k, downloadLocation = paste(getwd(), "sumstat_prot", sep = "/")) # Downloading the summary statistics for the protein of interest
  syn_code_path_list[[l]] <- syn_code$path
  syn_code_cache_list[[l]] <- syn_code$cacheDir
  l <- l+1
}

syn_code_cache <- getwd()

chr_list <- list('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')

###### MR ANALYSIS LOOP ######
j <- 1
for (j in 1:length(exposure_list)) {
  print(paste0("Starting on exposure: ", exposure_list[j]))
  
  mr_output <- NULL
  sensitivity_output <- NULL
  pleiotropy_output <- NULL
  
  instruments_final <- data.frame(1:5)
  colnames(instruments_final) = c("id")
  
  # Loading in the exposure summary statistics
  exposure_sumstats <- fread(exposure_file$path[j], header = TRUE, stringsAsFactors = FALSE)
  exposure_sumstats$CHR <- ifelse(exposure_sumstats$CHR == 23, "X", exposure_sumstats$CHR)
  
  # Processing adipose sum stats to filter by P value for significance
  exposure_sumstats <- exposure_sumstats[exposure_sumstats$P_BOLT_LMM < 5e-8]
  
  l <- 1
  for (i in sumstats_info$Code){
    
    print(paste0("Starting on outcomes: ", sumstats_info[sumstats_info$Code == i, ]$Assay[1]))
    log_message(paste0("Starting on outcome: ", sumstats_info[sumstats_info$Code == i, ]$Assay[1]))
    
    # Retrieving summary stats and reformatting
    untar(paste(syn_code_path_list[[l]]), list = F, exdir = paste(syn_code_cache))
    
    # Extracting all chromosome GWAS data for proteints
    chrom <- NULL
    for (m in chr_list){
      chrom_add <- fread(paste0(syn_code_cache, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code == i, ]$Docname[1]), "/", 
                                "discovery_chr", m, "_", sumstats_info[sumstats_info$Code == i, ]$UKBPPP_ProteinID[1], 
                                ":", sumstats_info[sumstats_info$Code == i, ]$Panel[1], ".gz"))
      chrom <- rbind(chrom, chrom_add)
      
    }
    
    l <- l+1
  
    chrom$P <- 10^-chrom$LOG10P
    chrom$CHROM <- ifelse(chrom$CHROM == 23, "X", chrom$CHROM) # Transferring 23rd chromosome to "X" for consistency
    chrom <- unique(chrom)
    
    if (is.null(chrom) || nrow(chrom) == 0) {
      print(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "chrom null"))
      log_message(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "chrom null"))
    } else{
      exposure_overlap <- exposure_sumstats[exposure_sumstats$GENPOS %in% chrom$GENPOS,] # Only selecting the variants that are overlapping
      exposure_overlap <- unique(exposure_overlap)
      
      if (is.null(exposure_overlap) || nrow(exposure_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "exposure overlap null"))
        log_message(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "exposure overlap null"))
      } else{
        exposure_rsid <- exposure_overlap[,c("CHR", "GENPOS", "SNP")] # Creating internal chr:pos -> rsid map to use later after harmonization but before clumping
        exposure_rsid <- exposure_rsid %>% mutate(SNP = strsplit(as.character(SNP), ",")) %>% unnest(SNP)
        exposure_overlap$phen <- paste(exposure_list[j])
        exposure_overlap$ID <- paste(exposure_overlap$CHR, exposure_overlap$GENPOS, exposure_overlap$ALLELE1, exposure_overlap$ALLELE0, sep = ":") # Harmonization using these unique id values to avoid prematurely filtering out duplicated SNPs
        exposure_overlap <- data.frame(exposure_overlap)
        exposure_overlap <- format_data(exposure_overlap, type=  "exposure", 
                                       phenotype_col = "phen", 
                                       snp_col = "ID", 
                                       beta_col = "BETA", 
                                       se_col = "SE", 
                                       eaf_col = "A1FREQ",
                                       effect_allele_col = "ALLELE1", 
                                       other_allele_col = "ALLELE0", 
                                       pval_col = "P_BOLT_LMM", 
                                       chr_col = "CHR", 
                                       pos_col = "GENPOS")
        
        # Finding overlapping chromosome and position in the outcome table
        chrom_overlap <- chrom[chrom$GENPOS %in% exposure_overlap$pos.exposure,] # Again, we just take the overlapping variants (now in the other direction)
        chrom_overlap_2 <- chrom_overlap # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
        chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
        chrom_overlap_2$A1FREQ <- 1 - chrom_overlap_2$A1FREQ
        colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <- c("ALLELE1", "ALLELE0")
        chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
        chrom_overlap$ID <- paste(chrom_overlap$CHROM, chrom_overlap$GENPOS, chrom_overlap$ALLELE1, chrom_overlap$ALLELE0, sep = ":") # Creating the same unique identifier as above for harmonization
        chrom_overlap$phen <- sumstats_info[sumstats_info$Code == i, ]$Assay[1]
        chrom_overlap <- data.frame(chrom_overlap)
        chrom_overlap <- format_data(chrom_overlap, type = "outcome", 
                                     phenotype_col = "phen", 
                                     snp_col = "ID", 
                                     beta_col = "BETA", 
                                     se_col = "SE", 
                                     eaf_col = "A1FREQ",
                                     effect_allele_col = "ALLELE1", 
                                     other_allele_col = "ALLELE0", 
                                     pval_col = "LOG10P", 
                                     chr_col = "CHROM", 
                                     samplesize_col = "N", 
                                     pos_col = "GENPOS", 
                                     log_pval = T)
        rm(chrom_overlap_2, chrom)
        
        # Harmonization
        dat <- harmonise_data(exposure_dat=exposure_overlap, outcome_dat=chrom_overlap) # This is where the matching happens
        
        # Reconciling custom ID and the Rsid values using GENPOS because all on same chromosome
        if (sumstats_info[sumstats_info$Code == i, ]$chr[1] == "X"){ # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
          dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x = "pos.outcome", by.y = "V2", all.x = T) 
          colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("CHR", "SNP")
        } else {
          dat <- merge(dat, exposure_rsid, by.x = "pos.outcome", by.y = "GENPOS", all.x = T) # Duplicates some rows since the chromosome and position are the same with different alleles
        }
        
        # Reformatting and organizing the data before clumping
        colnames(dat)[colnames(dat) %in% c("SNP.x", "SNP.y")] <- c("pos_id", "SNP") # We make sure that our reference column (which should have the name "SNP") is the RSID column
        dat <- dat[order(dat$pval.exposure), ] # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
        dat <- dat[!duplicated(dat$SNP), ] # Ordering the table by pval first ensures that the more significant SNP is kept
        dat <- dat[!duplicated(dat$pos_id), ]
        
        # Clumping
        # dat = clump_data(dat, clump_p1 = 5e-6, clump_p2 = 5e-6, clump_r2 = 0.1, clump_kb = 200)  # Now with the selection of the cis protein coding segment maybe need to make the window much smaller
        colnames(dat)[colnames(dat) %in% c("id.outcome", "id.exposure")] <- c("outcome_alt", "exposure_alt")
        dat$id.outcome <- dat$outcome
        dat$id.exposure <- dat$exposure
        
        
        # 2 clumping methods for MR-RAPS and MR with correlation matrix
        clump <- ld_clump(dplyr::tibble(rsid = dat$SNP, pval = dat$pval.exposure, id = dat$id.exposure), # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                          plink_bin = genetics.binaRies::get_plink_binary(), clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                          bfile = "/Users/michael.tian/Desktop/Natarajan_Lab/Tools/g1000_eur")
        dat2 <- dat[dat$SNP %in% clump$rsid, ] # dat2 for MR-Raps
        
        clump2 <- ld_clump(dplyr::tibble(rsid = dat$SNP, pval = dat$pval.exposure, id = dat$id.exposure), # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                           plink_bin = genetics.binaRies::get_plink_binary(), clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                           bfile = "/Users/michael.tian/Desktop/Natarajan_Lab/Tools/g1000_eur")
        dat <- dat[dat$SNP %in% clump2$rsid, ]
        
        rm(exposure_overlap, exposure_rsid, clump, clump2)
        
        # Steiger Filtering
        dat$samplesize.outcome <- chrom_overlap$samplesize.outcome[1]
        dat$samplesize.exposure <- exposure_file$samplesize[j]
        dat <- steiger_filtering(dat)
        dat <- dat[which(dat$steiger_dir == TRUE & dat$steiger_pval < 0.05), ]
        
        dat2$samplesize.outcome <- chrom_overlap$samplesize.outcome[1]
        dat2$samplesize.exposure <- exposure_file$samplesize[j]
        dat2 <- steiger_filtering(dat2)
        dat2 <- dat2[which(dat2$steiger_dir == TRUE & dat2$steiger_pval < 0.05), ]
        
        # MR Raps method
        res2 <- mr.raps(dat2$beta.exposure[dat2$mr_keep == T], dat2$beta.outcome[dat2$mr_keep == T], dat2$se.exposure[dat2$mr_keep == T], dat2$se.outcome[dat2$mr_keep == T])
        
        res2 <- data.frame(res2)
        res2$id.outcome <- sumstats_info[sumstats_info$Code == i, ]$Assay[1]
        res2$id.exposure <- exposure_list[j]
        res2$method <- "MR Raps"
        res2$nsnp <- nrow(dat2)
        colnames(res2)[colnames(res2) %in% c("beta.hat", "beta.se", "beta.p.value")] <- c("b", "se", "pval")
        res2 <- res2[c("id.exposure", "id.outcome", "nsnp", "method", "b", "se", "pval")]
        res2 <- generate_odds_ratios(res2)
        
        if (nrow(dat[dat$mr_keep,]) == 0) {
          print(paste0("Skipping ", sumstats_info[sumstats_info$Code= = i, ]$Assay[1], "mr_keep empty"))
          log_message(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "mr_keep empty"))
        } else {

          # Art's MR method accounting for between variant correlation with a lenient r2
          if (nrow(dat) == 1) { # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
            results_mr <- mr(dat, method_list = c("mr_wald_ratio"))
            res <- data.frame(exp = paste(exposure_list[j]), outc = sumstats_info[sumstats_info$Code == i, ]$Assay[1], nsnp = results_mr$nsnp, method = results_mr$method, b = results_mr$b,
                              se = results_mr$se, pval = results_mr$pval)
          } else if (nrow(dat) == 2) { # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
            ld <- ld_matrix(dat$SNP, bfile = "/Users/michael.tian/Desktop/Natarajan_Lab/Tools/g1000_eur", plink_bin = genetics.binaRies::get_plink_binary())
            rownames(ld) <- gsub("\\_.*", "", rownames(ld))
            dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
            dat3 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                     by = dat$beta.outcome, byse = dat$se.outcome,
                                                     correlation = ld)
            print(sumstats_info[sumstats_info$Code == i, ]$Assay[1])
            output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat3, correl = TRUE)
            res <- data.frame(exp = paste(exposure_list[j]), outc = sumstats_info[sumstats_info$Code == i, ]$Assay[1], nsnp = output_mr_ivw_corr@SNPs,
                              method = "Inverse variance weighted (correlation inc)", b = output_mr_ivw_corr@Estimate,
                              se = output_mr_ivw_corr@StdError, pval = output_mr_ivw_corr@Pvalue)
          } else { # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
            ld <- ld_matrix(dat$SNP, bfile = "/Users/michael.tian/Desktop/Natarajan_Lab/Tools/g1000_eur", plink_bin = genetics.binaRies::get_plink_binary())
            rownames(ld) <- gsub("\\_.*", "", rownames(ld))
            
            dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep = "_") %in% colnames(ld)), ]
            colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
            dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
            dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
            dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
            dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
            dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep = "_") %in% colnames(ld)), ]
            dat <- rbind(dat, dat_dupl)
            dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
            dat3 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                     by = dat$beta.outcome, byse = dat$se.outcome,
                                                     correlation = ld)
            
            print(sumstats_info[sumstats_info$Code==i,]$Assay[1])

            # IVW
            # Use tryCatch to handle errors in mr_ivw
            result <- tryCatch({
              output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat3, correl = TRUE)
              list(success = TRUE, result = output_mr_ivw_corr)
            }, error = function(e) {
              message("\n Error in mr_ivw: ", e$message)
              list(success = FALSE, result = NULL)
            })
    
            if (result$success) {
              results1 <- data.frame(exp = paste(exposure_list[j]), outc = sumstats_info[sumstats_info$Code == i,]$Assay[1], nsnp = result$result@SNPs,
                                     method = "Inverse variance weighted (correlation inc)", b = result$result@Estimate,
                                     se = result$result@StdError, pval = result$result@Pvalue)
            } else {
              # Handle the case when mr_ivw fails, e.g., set default values or skip processing
              cat(paste0("IVW failed for (IV >= 2) ", exposure_list[j], " ~ ", sumstats_info[which(sumstats_info$Code == i)[1], "Assay"]))
              log_message(paste0("IVW failed for (IV >= 2) ", exposure_list[j], " ~ ", sumstats_info[which(sumstats_info$Code == i)[1], "Assay"]))
              next
            }
            
            # MR-EGGER
            result2 <- tryCatch({
              output_mr_egger_corr <- MendelianRandomization::mr_egger(dat3, correl = TRUE)
              list(success = TRUE, result2 = output_mr_egger_corr)
            }, error = function(e) {
              message("\n Error in mr_egger: ", e$message)
              list(success = FALSE, result = NULL)
            })
            
            if (result2$success) {
              results2 <- data.frame(exp = paste(exposure_list[j]), outc = sumstats_info[sumstats_info$Code == i, ]$Assay[1], nsnp = output_mr_egger_corr@SNPs,
                                     method = "Egger (correlation inc)", b = output_mr_egger_corr@Estimate,
                                     se = output_mr_egger_corr@StdError.Est, pval = output_mr_egger_corr@Pvalue.Est)
              results3 <- data.frame(exp = paste(exposure_list[j]), outc = sumstats_info[sumstats_info$Code == i, ]$Assay[1], nsnp = output_mr_egger_corr@SNPs,
                                     method = "Egger intercept (correlation inc)", b = output_mr_egger_corr@Intercept,
                                     se = output_mr_egger_corr@StdError.Int, pval = output_mr_egger_corr@Pvalue.Int)
            } else {
              # Handle the case when mr_egger fails, e.g., set default values or skip processing
              cat(paste0("Egger failed for (IV >= 2) ", exposure_list[j], " ~ ", sumstats_info[which(sumstats_info$Code == i)[1], "Assay"]))
              log_message(paste0("Egger failed for (IV >= 2) ", exposure_list[j], " ~ ", sumstats_info[which(sumstats_info$Code == i)[1], "Assay"]))
              next
            }
            
            
            res <- rbind(results1, results2, results3)
            rm(results1, results2, results3)
          }
          res = generate_odds_ratios(res)
          colnames(res)[colnames(res) %in% c("exp", "outc")] <- c("id.exposure", "id.outcome")
          
          
          if (is.null(res) || nrow(res) == 0) {
            print(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "res empty"))
            log_message(paste0("Skipping ", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "res empty"))
          } else {
            # Merge onto all the MR results
            mr_output <- rbind(mr_output, res) 
            mr_output <- rbind(mr_output, res2) 
            
            # Sensitivity analysis
            sens <- mr_heterogeneity(dat)
            sensitivity_output <- rbind(sensitivity_output, sens) # all the sensitivity results
            
            sens2 <- mr_heterogeneity(dat2)
            sensitivity_output <- rbind(sensitivity_output, sens2) # all the sensitivity results
            
            # Horizontal pleiotropy test
            pleio <- mr_pleiotropy_test(dat)
            pleiotropy_output <- rbind(pleiotropy_output, pleio) # all the pleiotropy results
            
            pleio2 <- mr_pleiotropy_test(dat2)
            pleiotropy_output <- rbind(pleiotropy_output, pleio2) # all the pleiotropy results
            
            p1 <- mr_scatter_plot(res, dat)
            ggsave(p1[[1]], file = paste0(path, "/results/MR/", run, exposure_list[j], "/scatter_plot/", sumstats_info[sumstats_info$Code == i, ]$Assay[1], "_", exposure_list[j], ".pdf"), width = 7, height = 7)
          
            # Adding instruments to table
            snp_append <- unique(dat$SNP)
            instruments_new <- data.frame(1:length(snp_append), snp_append)
            colnames(instruments_new) <- c("id", paste0(dat$outcome[1]))
            instruments_final <- merge(instruments_final, instruments_new, by = "id", all = TRUE)
          }
        }
      }
    }
  }
  
  
  # Saving all of these outputs per outcome
  mr_output <- as.data.frame(mr_output)
  sensitivity_output <- as.data.frame(sensitivity_output)
  pleiotropy_output <- as.data.frame(pleiotropy_output)
  
  all_outcome_results <- rbind(all_outcome_results, mr_output)
  
  # Writing out files: 1 per outcome
  fwrite(mr_output, paste0(path, "/results/MR/", run, exposure_list[j], "/mr_results_", exposure_list[j], "_", date, ".csv")) 
  fwrite(sensitivity_output, paste0(path, "/results/MR/", run, exposure_list[j], "/mr_sensitivity_", exposure_list[j], "_", date, ".csv")) 
  fwrite(pleiotropy_output, paste0(path, "/results/MR/", run, exposure_list[j], "/mr_pleiotropy_", exposure_list[j], "_", date, ".csv"))
  
  instruments_final <- as.data.frame(instruments_final)
  fwrite(instruments_final, paste0(path, "/results/MR/", run, "genetic_instruments_", exposure_list[j], "_", date, ".csv"))
}

# Saving the results from all outcomes altogether
fwrite(all_outcome_results, paste0(path, "/results/MR/", run, "all_outcomes_mr_results_", date, ".csv")) 

close(log_connection)
