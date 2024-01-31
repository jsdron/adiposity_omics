#### 0 - LIBRARIES ####

# only need to do these once. comment them out and save after they have been successfully run
install.packages("remotes")
library(remotes) 
install_github("MRCIEU/TwoSampleMR") # only do this once. after, you can just load the library

# run these every time
library(TwoSampleMR)
library(data.table)
library(tidyr)
library(dplyr)


#### 1 - LOAD SUMMARY STATISTICS META-DATA ####

# adiposity outcomes
outcome_ss <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/adiposity_sumstats.txt') # JD created this file

# the exposures are the various analytes (metabolomic or proteomic)
# metabolites
exposure_ss <- fread('') # JD created this file

#### 2 - RUN MR ####
# Run MR, sensitivity analysis, and horizontal pleiotropy test for each outcome

all_outcome_results = NULL

for (j in 1:length(outcome_ss)) { # will loop through one outcome at a time 
 
  print("Loading summary statistics...")
  
  print(paste0("Formatting outcome:", outcome_ss$variable_name[j] ))

  outcome_path <- outcome_ss$path[j]
  outcome <- fread(outcome_path)
  outcome <- subset(outcome, outcome$INFO >= 0.7)
  outcome$phenotype <- outcome_ss$variable_name[j] 
  
  outcome <- format_data(outcome, type = "outcome",
                        snp_col = "SNP", # combo of rs367896724 and 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C
                        beta_col = "BETA", # effect size from BOLT-LMM approximation to infinitesimal mixed model 
                        se_col = "SE",
                        eaf_col = "A1FREQ",
                        effect_allele_col = "ALLELE1",
                        other_allele_col = "ALLELE0",
                        pval_col = "P_BOLT_LMM", # non-infinitesimal mixed model association test p-value
                        phenotype_col = "phenotype")
  
  mr_output = NULL
  sensitivity_output = NULL
  pleiotropy_output = NULL
  all_dat = NULL
  
  for (i in 1:length(exposure_ss)) { # will loop through one analyte at a time
   
    print(paste0(i, ": ", exposure_ss$NAME[i], " for ", j, ": ", outcome_ss$variable_name[j]))
    
    ### Forward MR: Exposure -> Outcome ###
    # Grab the significant exposure SNPs and then filter against the outcome summary statistics, whether the outcome SNPs are sig or not
    exposure <- read_exposure_data(filename,
      sep = " ",
      phenotype_col = exposure_ss$NAME[i],
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      units_col = "units",
      ncase_col = "ncase",
      ncontrol_col = "ncontrol",
      samplesize_col = "samplesize",
      gene_col = "gene",
      id_col = "id",
      min_pval = 5e-8, # genome-wide significance
      log_pval = FALSE,
      chr_col = "chr",
      pos_col = "pos"
    )
    
    dat <- harmonise_data(exposure, outcome)
    dat <- clump_data(dat) # defaults here: https://mrcieu.github.io/TwoSampleMR/reference/clump_data.html
    all_dat <- rbind(all_dat, dat) # all the data pre-MR
    
    res = mr(dat) 
    ### res = generate_odds_ratios(res)
    res$exposure <- res$id.exposure
    mr_output = rbind(mr_output, res) # all the MR results
    
    ### Sensitivity analysis ###
    sens <- mr_heterogeneity(dat)
    sensitivity_output = rbind(sensitivity_output, sens) # all the sensitivity results
    
    ### Horizontal pleiotropy test ###
    pleio <- mr_pleiotropy_test(dat)
    pleiotropy_output = rbind(pleiotropy_output, pleio) # all the pleiotropy results
    
    p1 <- mr_scatter_plot(res, dat)
    ggsave(p1[[1]], file = paste0(path, "/results/", run, outcome_list[j],"/scatter_plot/",outcome_list[j],"_",exposure_names[i],".pdf"), width = 7, height = 7)
    
  }
  mr_output = as.data.frame(mr_output)
  sensitivity_output = as.data.frame(sensitivity_output)
  pleiotropy_output = as.data.frame(pleiotropy_output)
  
  # Writing out files: 1 per outcome
  fwrite(mr_output, paste0(path, "/results/", run, outcome_list[j], "/mr_results_",outcome_list[j],"_",date,".csv")) 
  fwrite(sensitivity_output, paste0(path, "/results/", run, outcome_list[j], "/mr_sensitivity_",outcome_list[j],"_",date,".csv")) 
  fwrite(pleiotropy_output, paste0(path, "/results/", run, outcome_list[j], "/mr_pleiotropy_",outcome_list[j],"_",date,".csv")) 
  
  all_outcome_results <- rbind(all_outcome_results, mr_output)
  
  outcome_sumstats = NULL
  outcome = NULL
}

fwrite(all_outcome_results, paste0(path, "/results/", run, "all_outcomes_mr_results_",date,".csv")) 

# Extract the genetic instruments
exposure = extract_instruments(outcomes = exposure_list[1])
instruments_final = data.frame(1:nrow(exposure), exposure$SNP)
colnames(instruments_final) = c("id", paste0(exposure_list[1]))

for (i in 2:length(exposure_list)) {
  ## for (i in 2:5) {
  print(i)
  exposure = extract_instruments(outcomes = exposure_list[i])
  exposure_new = data.frame(1:nrow(exposure), exposure$SNP)
  colnames(exposure_new) = c("id", paste0(exposure_list[i]))
  instruments_final = merge(instruments_final, exposure_new, by = "id", all = TRUE)
}

instruments_final = as.data.frame(instruments_final)
fwrite(instruments_final, paste0(path, "/results/", run, "genetic_instruments_for_exposures_",date,".csv")) 

