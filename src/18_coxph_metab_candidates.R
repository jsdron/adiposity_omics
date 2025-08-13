#########################################
# Script: coxph_metab_candidates.R
# Description: Fits Cox proportional hazards models for metabolite candidates 
#              associated with ASAT, VAT, and GFAT, and generates forest plots 
#              for their associations with coronary artery disease, type 2 diabetes, and chronic kidney disease.
# Key Outputs:
#   - Figure (forest plots for metabolite associations with CAD, T2D, CKD)
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
metab_exposures <- c("HDL_CE", "M_HDL_CE", "M_HDL_C", "HDL_L", "HDL_PL", "ApoA1", "M_HDL_FC","M_HDL_P",
               "M_HDL_L", "HDL_P", "M_HDL_PL")  # GFAT 
event_names <- c('t2d', 'cad', 'ckd') 
covar <- c("enroll_age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")  
covariates <- paste(covar, collapse = " + ")  

threshold <- 0.05/length(metab_exposures)

###### PREPARE DATAFRAME ######
metab <- data.frame(fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb676832/ukb676832.ukbnmr.tsv.gz',
                          select=c("eid", "HDL_CE", "M_HDL_CE", "M_HDL_C", "HDL_L", "HDL_PL", "ApoA1", "M_HDL_FC","M_HDL_P",
                                   "M_HDL_L", "HDL_P", "M_HDL_PL")))

diseases <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/bri/input/BRI_epi.2025-02-05.tsv", 
                  select=c("eid", "sex", "enroll_age", "BMI" , "knn",  
                           "incd_t2d", "prev_t2d", "censor_age_t2d",
                           "incd_ckd", "prev_ckd", "censor_age_ckd",
                           "incd_cad", "prev_cad", "censor_age_cad",
                           "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

main <- merge(diseases, metab, by="eid") # 

withdrawn <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/withdrawn_samples/w7089_2023-08-21.csv')
colnames(withdrawn)[1] <- "eid"
withdrawn_vector <- c(withdrawn$eid)
main <- main[!(main$eid %in% withdrawn_vector),]
rm(withdrawn, withdrawn_vector)

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

for (i in (which(colnames(main) == "HDL_CE"):which(colnames(main) == "M_HDL_PL"))) {
  main[,colnames(main)[i]] <- scale(as.numeric(as.matrix(main[,..i])))
}

###### COXPH MODEL ######
results_df <- data.frame()  

for (i in 1:length(event_names)){
  
  disease_name <- event_names[i]  
  incd_event <- paste0("incd_", event_names[i])  
  
  # Loop through each exposure and run a separate model for each
  for (exposure in metab_exposures) {
    
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
      mutate(p.value = summary(model)$coefficients[1, "Pr(>|z|)"]) %>%
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

results_df$depot <- "GFAT"

results_df <- results_df %>% distinct()

# Save results
write.table(results_df, file="/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/coxph/coxph_metab.2025-05-07.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)

sig <- subset(results_df, significance=="Significant")

###### FOREST PLOTS ######
results_metab <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/save_metabolites.tsv", data.table = F)
metab_exposures <- c("HDL_CE", "M_HDL_CE", "M_HDL_C", "HDL_L", "HDL_PL", "ApoA1", "M_HDL_FC","M_HDL_P",
                     "M_HDL_L", "HDL_P", "M_HDL_PL")

sum(duplicated(results_metab)) # 0 duplicates

threshold <- 0.05/length(metab_exposures)
results_metab$significance <- ifelse(results_metab$pval < threshold, "Significant", "Not Significant")
results_metab$depot <- "GFAT"

results_metab <- results_metab %>% distinct()

results_metab$omic <- "Metabolomic"

results_metab$N_label <- paste0(results_metab$n_incidence, " / ", results_metab$n)

# Subsetting and creating meta-analysis objects
results_cad_metab <- subset(results_metab, results_metab$outcome=="Coronary_Artery_Disease")
results_t2d_metab <- subset(results_metab, results_metab$outcome=="Diabetes_Type_2")
results_ckd_metab <- subset(results_metab, results_metab$outcome=="Chronic_kidney_disease")

# CAD ----
meta_analysis_cad_metab <- metagen(
  TE = log(HR), seTE = se,
  studlab = exposures, byvar = depot, sm = "HR", 
  lower = log(lci), upper = log(uci),
  n.c = n, n.e = n_incidence, 
  data = results_cad_metab,
  pval = pval)

meta_analysis_cad_metab$N_label <- results_cad_metab$N_label
meta_analysis_cad_metab$my_pval <- sprintf("%.2e", meta_analysis_cad_metab$pval)

pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/jd_version/Figure5_metab_cad.pdf",
    width = 2.5,              # inches; adjust as needed
    height = 2.5,  # simple heuristic
    pointsize = 7,            # <-- sets base font to 7 pt
    family = "Arial",
    useDingbats = FALSE)      # safer glyph embedding for journals

forest(
  x = meta_analysis_cad_metab,
  common = FALSE, random = FALSE, overall = FALSE,
  col.square=c("#E1D2B2"), col.square.lines=c("#E1D2B2"),
  sortvar = TE,
  print.subgroup.hetstat = FALSE,
  print.overall.hetstat = FALSE,
  fontsize=7,
  weight.study = "same",
  colgap = unit(4, "mm"),                 # narrower column gaps
  colgap.forest.left = "4mm",              # less gap before forest plot
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
meta_analysis_t2d_metab <- metagen(
  TE = log(HR), seTE = se,
  studlab = exposures, byvar = depot, sm = "HR", 
  lower = log(lci), upper = log(uci),
  n.c = n, n.e = n_incidence, 
  data = results_t2d_metab,
  pval = pval)

meta_analysis_t2d_metab$N_label <- results_t2d_metab$N_label
meta_analysis_t2d_metab$my_pval <- sprintf("%.2e", meta_analysis_t2d_metab$pval) 

pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/jd_version/Figure5_metab_t2d.pdf",
    width = 2.5,              # inches; adjust as needed
    height = 2.5,  # simple heuristic
    pointsize = 7,            # <-- sets base font to 7 pt
    family = "Arial",
    useDingbats = FALSE)      # safer glyph embedding for journals

forest(
  x = meta_analysis_t2d_metab,
  common = FALSE, random = FALSE, overall = FALSE,
  col.square=c("#E1D2B2"), col.square.lines=c("#E1D2B2"),
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
meta_analysis_ckd_metab <- metagen(
  TE = log(HR), seTE = se,
  studlab = exposures, byvar = depot, sm = "HR", 
  lower = log(lci), upper = log(uci),
  n.c = n, n.e = n_incidence, 
  data = results_ckd_metab,
  pval = pval)

meta_analysis_ckd_metab$N_label <- results_ckd_metab$N_label
meta_analysis_ckd_metab$my_pval <- sprintf("%.2e", meta_analysis_ckd_metab$pval) 

pdf("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/linke_cox_model/jd_version/Figure5_metab_ckd.pdf",
    width = 2.5,              # inches; adjust as needed
    height = 2.5,  # simple heuristic
    pointsize = 7,            # <-- sets base font to 7 pt
    family = "Arial",
    useDingbats = FALSE)      # safer glyph embedding for journals

forest(
  x = meta_analysis_ckd_metab,
  common = FALSE, random = FALSE, overall = FALSE,
  col.square=c("#E1D2B2"), col.square.lines=c("#E1D2B2"),
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

