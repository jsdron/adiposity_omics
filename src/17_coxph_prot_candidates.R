#########################################
# Script: coxph_prot_candidates.R
# Description: Fits Cox proportional hazards models for protein candidates 
#              associated with ASAT, VAT, and GFAT, and generates forest plots 
#              for their associations with coronary artery disease, type 2 diabetes, and chronic kidney disease.
# Key Outputs:
#   - Figure (forest plots for protein associations with CAD, T2D, and CKD)
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(broom)
library(meta)
library(cmprsk)
library(survival)
library(survminer)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(extrafont)

font_import()  # This imports all available fonts to R
loadfonts(device = "pdf")  # Load fonts for PDF output

###### SET VARIABLES ######
exposures <- c("THBS2", # ASAT (FWD)
               "ABL1", "CCL17", "MSR1", "CFB", "KHK", "TSPAN8", "SELPLG", "ANXA2", "ASAH2", "ADAMTSL5", "SHBG", "ITGB6", # VAT (FWD)
               "SHBG", "LPL", "IL2RA", "TREH", "FCAMR", "NFASC") # GFAT (FWD)
event_names <- c('t2d', 'cad', 'ckd') 
covar <- c("enroll_age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")  
covariates <- paste(covar, collapse = " + ")  

threshold <- 0.05/length(exposures)

###### PREPARE DATAFRAME ######
prot <- data.frame(fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb672544/olink_data3k.txt.gz')) 
names_prot <- data.frame(fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb672544/coding143.tsv')) # Load linker file to assign protein names
names_prot <- separate(names_prot, meaning, into = c("short_name", "long_name"), sep = ";")

# Select only those with baseline protein measurements
prot <- prot[prot$ins_index == 0, !(names(prot) %in% "ins_index")]

# Make dataframe wide
prot <- as.data.frame(pivot_wider(prot, names_from = protein_id, values_from = result))
names_prot$coding <- as.character(names_prot$coding)
setnames(prot, old = names_prot$coding, new = names_prot$short_name, skip_absent = TRUE)

# Only keep the relevant columns
prot <- prot[, colnames(prot) %in% c("eid", exposures)]

diseases <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/bri/input/BRI_epi.2025-02-05.tsv", 
                  select=c("eid", "sex", "enroll_age", "BMI" , "knn",  
                           "incd_t2d", "prev_t2d", "censor_age_t2d",
                           "incd_ckd", "prev_ckd", "censor_age_ckd",
                           "incd_cad", "prev_cad", "censor_age_cad", 
                           "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

main <- merge(prot, diseases, by = "eid") # 
rm(prot, diseases)

withdrawn <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/withdrawn_samples/w7089_2023-08-21.csv')
colnames(withdrawn)[1] <- "eid"
withdrawn_vector <- c(withdrawn$eid)
main <- main[!(main$eid %in% withdrawn_vector),]
rm(withdrawn, withdrawn_vector)

# Restrict to white only
#main <- subset(main, knn == "EUR")

main <- main %>%  
  filter(if_all(all_of(covar), ~ !is.na(.))) 

no_prev_cad <- subset(main, prev_cad != 1) 
no_prev_t2d <- subset(main, main$prev_t2d != 1) 
no_prev_ckd <- subset(main, main$prev_ckd != 1) 
main <- subset(main, prev_cad != 1 & prev_t2d != 1 & prev_t2d != 1) 

table(main$sex) 
mean(main$enroll_age) 
sd(main$enroll_age) 

table(main$incd_cad) 
median(main$censor_age_cad - main$enroll_age) 
IQR(main$censor_age_cad - main$enroll_age) 

table(main$incd_t2d)
median(main$censor_age_t2d - main$enroll_age) 
IQR(main$censor_age_t2d - main$enroll_age) 

main <- data.frame(main)

for (i in (which(colnames(main) == "FCAMR"):which(colnames(main) == "ANXA2"))) {
  main[,colnames(main)[i]] <- scale(as.numeric(as.matrix(main[,i])))
}

###### COXPH MODEL ######
results_df <- data.frame()  

for (i in 1:length(event_names)){
  
  disease_name <- event_names[i]  
  incd_event <- paste0("incd_", event_names[i])  

  # Loop through each exposure and run a separate model for each
  for (exposure in exposures) {
    
    # Remove people with missing important variables  
    no_prev <- main %>%  
      filter(if_all(all_of(exposure), ~ !is.na(.))) %>%  
      filter(if_all(all_of(covar), ~ !is.na(.)))  
    
    # Using age of enrollment  
    no_prev$time <- no_prev[[paste0("censor_age_", disease_name)]] 
    
    
    # Create formula for each exposure
    exposure_formula <- as.formula(paste0("Surv(time, ", incd_event, ") ~ ", exposure, " + ", covariates))  
    
    # Run the Cox model for the current exposure
    model <- coxph(exposure_formula, data = no_prev)  
    
    # Extract and tidy the model results
    tidy_results <- tidy(model, conf.int = TRUE) %>%  
      mutate(term = exposure, 
             HR = exp(estimate), 
             lower_CI = exp(conf.low), 
             upper_CI = exp(conf.high)) %>%  
      slice(1)  # Keep only the first row with the effect for the current exposure  
    
    # Add outcome and sample size details
    tidy_results <- tidy_results %>%  
      mutate(outcome = disease_name, 
             N = length(no_prev$eid), 
             N.e = nrow(subset(no_prev, no_prev[[incd_event]] == 1)))
    
    tidy_results$mean_fu <- paste0(round(mean(no_prev$time), 2), " (", round(sd(no_prev$time), 2), ")")
    tidy_results$male <- sum(no_prev$sex)
    tidy_results$age <- paste0(round(mean(no_prev$enroll_age), 2), " (", round(sd(no_prev$enroll_age), 2), ")")
    tidy_results$event_per <- paste0(tidy_results$N.e, " (", round((tidy_results$N.e/tidy_results$N)*100, 2), "%)")
    
    # Append the results to the final dataframe  
    results_df <- bind_rows(results_df, tidy_results)  
  }
}

# Remove duplicates
results_df <- results_df %>%
  group_by(term, outcome) %>%
  filter(N == max(N)) %>%
  ungroup()

results_df$N_label <- paste0(results_df$N.e, " / ", results_df$N)

results_df$significance <- ifelse(results_df$p.value < threshold, "Significant", "Not Significant")
results_df$outcome <- ifelse(results_df$outcome=="cad", "Coronary Artery Disease", 
                             ifelse(results_df$outcome=="ckd", "Chronic Kidney Disease","Type 2 Diabetes"))

results_df$depot <- ifelse(results_df$term %in% c("THBS2"), "ASAT",
                           ifelse(results_df$term %in% c("LPL", "IL2RA", "TREH", "FCAMR", "NFASC"), "GFAT", 
                                  ifelse(results_df$term %in% c("SHBG"), "GFAT and VAT", 
                                         #ifelse(results_df$term %in% c("ITIH3", "SHBG","INHBC"), "GFAT (both) & VAT (fwd)", 
                                                "VAT")))#)

results_df <- results_df %>% distinct()

# Save results
write.table(results_df, file = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/coxph/coxph_prot.2025-05-07.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

###### FOREST PLOTS ######
results_prot <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/save.tsv", data.table = F)

results_prot$depot <- case_when(
  results_prot$exposures %in% c("THBS2")~"ASAT", 
  results_prot$exposures %in% c("LPL", "IL2RA", "TREH", "FCAMR", "NFASC") ~ "GFAT",
  results_prot$exposures %in% c("SHBG") ~ "GFAT and VAT",
  results_prot$exposures %in% c("ABL1","CCL17","MSR1", "CFB", "KHK", "TSPAN8", "SELPLG", "ANXA2", "ASAH2", "ADAMTSL5", "ITGB6") ~ "VAT",
  TRUE ~ NA
)

sum(duplicated(results_prot))
dup <- results_prot[duplicated(results_prot) | duplicated(results_prot, fromLast = T),]

## removing duplicates
results_prot <- results_prot %>%
  filter(!duplicated(.))

threshold <- 0.05/length(unique(results_prot$exposures))
results_prot$significance <- ifelse(results_prot$pval < threshold, "Significant", "Not Significant")

results_prot$omic <- "Proteomic"

results_prot$N_label <- paste0(results_prot$n_incidence, " / ", results_prot$n)

# Subsetting and creating meta-analysis objects
results_cad_prot <- subset(results_prot, results_prot$outcome=="Coronary_Artery_Disease")
results_t2d_prot <- subset(results_prot, results_prot$outcome=="Diabetes_Type_2")
results_ckd_prot <- subset(results_prot, results_prot$outcome=="Chronic_kidney_disease")

# CAD ----
meta_analysis_cad_prot <- metagen(
  TE = log(HR), seTE = se,
  studlab = exposures, byvar = depot, sm = "HR", 
  lower = log(lci), upper = log(uci),
  n.c = n, n.e = n_incidence, 
  data = results_cad_prot,
  pval = pval)

meta_analysis_cad_prot$N_label <- results_cad_prot$N_label
meta_analysis_cad_prot$my_pval <- sprintf("%.2e", meta_analysis_cad_prot$pval) 

pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/jd_version/Figure5_prot_cad.pdf",
    width = 2.5,              # inches; adjust as needed
    height = 4,  # simple heuristic
    pointsize = 7,            # <-- sets base font to 7 pt
    family = "Arial",
    useDingbats = FALSE)      # safer glyph embedding for journals

forest(
  x = meta_analysis_cad_prot,
  common = FALSE, random = FALSE, overall = FALSE,
  col.square=c("#A68A6D"), col.square.lines=c("#A68A6D"),
  sortvar = TE,
  print.subgroup.hetstat = FALSE,
  print.overall.hetstat = FALSE,
  fontsize=7,
  weight.study = "same",
  colgap = unit(4, "mm"),                 # narrower column gaps
  colgap.forest.left = "5mm",              # less gap before forest plot
  plotwidth = unit(2, "cm"),             # shrink forest panel width
  big.mark = ",",
  subgroup = FALSE, print.subgroup.name = FALSE, col.subgroup = "black",
  smlab = "", smlab.pos = 0,
  xlab = "Hazard Ratio (95% CI)",
  leftcols = c("studlab"), #, "N_label"),
  leftlabs = c("Exposure"), #, "Number of events /\nTotal participants"),
  rightcols = c("my_pval"), # "effect", "ci", 
  rightlabs = c("P-value"), # "HR", "95% CI",
  cex = 0.8                                # slight downscale of text/columns
)

dev.off()

# T2D ----
meta_analysis_t2d_prot <- metagen(
  TE = log(HR), seTE = se,
  studlab = exposures, byvar = depot, sm = "HR", 
  lower = log(lci), upper = log(uci),
  n.c = n, n.e = n_incidence, 
  data = results_t2d_prot,
  pval = pval)

meta_analysis_t2d_prot$N_label <- results_t2d_prot$N_label
meta_analysis_t2d_prot$my_pval <- sprintf("%.2e", meta_analysis_t2d_prot$pval) 

pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/jd_version/Figure5_prot_t2d.pdf",
    width = 2.5,              # inches; adjust as needed
    height = 4,  # simple heuristic
    pointsize = 7,            # <-- sets base font to 7 pt
    family = "Arial",
    useDingbats = FALSE)      # safer glyph embedding for journals

forest(
  x = meta_analysis_t2d_prot,
  common = FALSE, random = FALSE, overall = FALSE,
  col.square=c("#A68A6D"), col.square.lines=c("#A68A6D"),
  sortvar = TE,
  print.subgroup.hetstat = FALSE,
  print.overall.hetstat = FALSE,
  fontsize=7,
  weight.study = "same",
  colgap = unit(4, "mm"),                 # narrower column gaps
  colgap.forest.left = "5mm",              # less gap before forest plot
  plotwidth = unit(2, "cm"),             # shrink forest panel width
  big.mark = ",",
  subgroup = FALSE, print.subgroup.name = FALSE, col.subgroup = "black",
  smlab = "", smlab.pos = 0,
  xlab = "Hazard Ratio (95% CI)",
  leftcols = c("studlab"), #, "N_label"),
  leftlabs = c("Exposure"), #, "Number of events /\nTotal participants"),
  rightcols = c("my_pval"), # "effect", "ci", 
  rightlabs = c("P-value"), # "HR", "95% CI",
  cex = 0.8                                # slight downscale of text/columns
)

dev.off()

# CKD ----
meta_analysis_ckd_prot <- metagen(
  TE = log(HR), seTE = se,
  studlab = exposures, byvar = depot, sm = "HR", 
  lower = log(lci), upper = log(uci),
  n.c = n, n.e = n_incidence, 
  data = results_ckd_prot,
  pval = pval)

meta_analysis_ckd_prot$N_label <- results_ckd_prot$N_label
meta_analysis_ckd_prot$my_pval <- sprintf("%.2e", meta_analysis_ckd_prot$pval) 

pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/jd_version/Figure5_prot_ckd.pdf",
    width = 2.5,              # inches; adjust as needed
    height = 4,  # simple heuristic
    pointsize = 7,            # <-- sets base font to 7 pt
    family = "Arial",
    useDingbats = FALSE)      # safer glyph embedding for journals

forest(
  x = meta_analysis_ckd_prot,
  common = FALSE, random = FALSE, overall = FALSE,
  col.square=c("#A68A6D"), col.square.lines=c("#A68A6D"),
  sortvar = TE,
  print.subgroup.hetstat = FALSE,
  print.overall.hetstat = FALSE,
  fontsize=7,
  weight.study = "same",
  colgap = unit(4, "mm"),                 # narrower column gaps
  colgap.forest.left = "5mm",              # less gap before forest plot
  plotwidth = unit(2, "cm"),             # shrink forest panel width
  big.mark = ",",
  subgroup = FALSE, print.subgroup.name = FALSE, col.subgroup = "black",
  smlab = "", smlab.pos = 0,
  xlab = "Hazard Ratio (95% CI)",
  leftcols = c("studlab"), #, "N_label"),
  leftlabs = c("Exposure"), #, "Number of events /\nTotal participants"),
  rightcols = c("my_pval"), # "effect", "ci", 
  rightlabs = c("P-value"), # "HR", "95% CI",
  cex = 0.8                                # slight downscale of text/columns
)

dev.off()
