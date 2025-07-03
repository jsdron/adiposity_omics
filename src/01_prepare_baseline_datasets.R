#########################################
# Script: prepare_baseline_datasets.R
# Description: Prepares and cleans baseline data with covariates, 
#              adiposity measures, proteomics and metabolomics data,
#              removes withdrawn samples and relatives, imputes and scales omics variables,
#              and saves cleaned datasets.
# Key Output: Datasets to use for the epidemiological analyses of the paper.
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyr)
library(dplyr)
library(impute)

###### PREPARE DATAFRAME(S) ######
# Load baseline covariates 
baseline <- fread("/Volumes/medpop_esp2/aniruddh/Lpa/UKBB_TG_BackgroundVariables_22020.txt",
                  select = c("eid", "enroll_age", "sex", "ethnicity", "eversmoker", "eversmoked", "currentsmoker", "DBP", "SBP", "fasting1", "anylipidmed0", "cholmed0", "statin0")
                  )

tmp <- fread("/Volumes/medpop_esp2/aschuerm/resources/ukb_df_base_full_20230905.tsv.gz",  
             select = c("id", "BMI", "race", "antihtnbase", "BP_Meds")
             )

# Merge and keep environment tidy
baseline <- merge(baseline, tmp, by.x = "eid", by.y = "id")
rm(tmp)

# Create fasting status and BMI category
baseline$fasting_status <- ifelse(baseline$fasting1 >= 8, 1, 0) # considered fasted (1) if fasting time is greater than or equal to 8 hours
baseline$BMI_class <- ifelse(baseline$BMI < 18.5, "underweight",
                             ifelse(baseline$BMI >= 18.5 & baseline$BMI < 25, "normal_weight",
                                    ifelse(baseline$BMI >= 25 & baseline$BMI < 30, "overweight", "obesity")
                                    )
                             )

# Recode the sex variable
baseline$sex <- ifelse(baseline$sex == 1, "M", 
                         ifelse(baseline$sex == 0, "F", NA))

# Rename columns
colnames(baseline) <- c(
  "eid", "enroll_age", "sex", "ethnicity", "eversmoker", "eversmoked", "currentsmoker",
  "DBP", "SBP", "fasting_time", "anylipidmed", "chol_med", "statin", 
  "BMI", "race", "HTN_med", "BP_med", "fasting_status", "BMI_class"
)

# Select and order columns
baseline <- baseline[, c(
  "eid", "enroll_age", "sex", "ethnicity", "race", "eversmoker", "eversmoked", "currentsmoker",
  "BMI", "BMI_class", "DBP", "SBP", "fasting_time", "fasting_status", 
  "anylipidmed", "chol_med", "statin", "HTN_med", "BP_med"
)]

# Load MRI-derived adiposity measurements
adi <- fread("/Volumes/medpop_esp2/sagrawal/shared_fat/data/fatPhenoCovar_wFMR.csv")
adi <- adi[, -c("sex")]

# Merge MRI data
baseline <- merge(baseline, adi, by = "eid", all.y = TRUE)
baseline$time_between <- baseline$age_instance2 - baseline$enroll_age
rm(adi)

# Remove withdrawn participants
withdrawn <- fread("/Volumes/medpop_esp2/projects/UK_Biobank/withdrawn_samples/w7089_2023-08-21.csv")
colnames(withdrawn)[1] <- "eid"
withdrawn_vector <- c(withdrawn$eid)
baseline <- baseline[!(baseline$eid %in% withdrawn_vector), ]
rm(withdrawn, withdrawn_vector)

# Keep participants who self-reported "white"
baseline <- subset(baseline, baseline$race=="white")

# Create backup dataframe
save <- baseline # so we can use later to see how many ppl don't have omics measurements

###### METABOLOMICS DATA ######
# Load QC'd NMR metabolomics data
qc_met <- fread("/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb676832/ukb676832.ukbnmr.tsv.gz")

# Identify duplicated IDs (participants with multiple visits)
double_measures <- qc_met$eid[duplicated(qc_met$eid)]

# Keep only baseline measurements (visit_index == 0)
qc_met <- subset(qc_met, qc_met$visit_index == 0)

# Subset columns for analytes
index <- c(
  which(colnames(qc_met) == "eid"),
  which(colnames(qc_met) == "Clinical_LDL_C"):which(colnames(qc_met) == "Omega_6_pct_PUFA")
  )
b <- qc_met[, ..index]

# Remove features and samples with >10% missingness
b <- b %>% select(where(~mean(is.na(.)) <= 0.1))
b <- b[rowMeans(is.na(b[, -1])) <= 0.1, ]

# Keep only participants passing QC
qc_met <- qc_met[qc_met$eid %in% b$eid, ]

# Annotate metabolomics presence in baseline
baseline$metab <- ifelse(baseline$eid %in% qc_met$eid, 1, 0)
baseline <- merge(baseline, qc_met, by = "eid", all.x = TRUE)
baseline$double_metab <- ifelse(baseline$eid %in% double_measures, 1, 0)

# Also update save object for comparison
save$metab <- ifelse(save$eid %in% qc_met$eid, 1, 0)

rm(qc_met, double_measures, index, b)

###### PROTEOMICS DATA ######
# Load long-format proteomics data
prot <- fread("/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb672544/olink_data3k.txt.gz")

# Load protein ID-to-name mapping
names_prot <- fread("/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb672544/coding143.tsv")
names_prot <- separate(names_prot, meaning, into = c("short_name", "long_name"), sep = ";")

# Restrict to baseline and pivot to wide format
prot <- prot[prot$ins_index == 0, -("ins_index")]
prot <- as.data.frame(pivot_wider(prot, names_from = protein_id, values_from = result))

# Map coding to short protein names
names_prot$coding <- as.character(names_prot$coding)
setnames(prot, old = names_prot$coding, new = names_prot$short_name, skip_absent = TRUE)

# Identify protein analyte range
first_prot <- colnames(prot)[2]
last_prot <- colnames(prot)[ncol(prot)]

# Annotate proteomics presence
baseline$prot <- ifelse(baseline$eid %in% prot$eid, 1, 0)
table(baseline$prot)

# Restrict to protein columns for missingness filtering
index <- c(
  which(colnames(prot) == "eid"),
  which(colnames(prot) == first_prot):which(colnames(prot) == last_prot)
  )
b <- prot[, index]

# Remove proteins with >20% missingness
c <- b %>% select(where(~mean(is.na(.)) <= 0.2))
removed_proteins <- setdiff(names(b), names(c))

# Restrict to individuals in baseline and drop removed proteins
c <- subset(c, c$eid %in% baseline$eid)
c <- c[, !names(c) %in% removed_proteins]

# Merge cleaned proteomics data with baseline
baseline$prot <- ifelse(baseline$eid %in% c$eid, 1, 0)
baseline <- merge(baseline, c, by = "eid", all.x = TRUE)

# Final sample restriction to participants with any omics
baseline <- subset(baseline, baseline$prot == 1 | baseline$metab == 1)

# Define omics group
baseline$omics <- ifelse(baseline$metab == 1 & baseline$prot == 1, "both", 
                         ifelse(baseline$metab == 1, "metab", "prot"))

# Update tracking dataframe
save$prot <- ifelse(save$eid %in% c$eid, 1, 0)
save$no_omics <- ifelse(save$metab == 0 & save$prot == 0, "none", "some")

# Check number of participants without omics
table(save$no_omics)

rm(prot, names_prot, index, b, c)

###### REMOVE RELATIVES ######
# Load relatedness data
related <- fread("/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb7089_rel_s488264.dat")

# Keep pairs with kinship ≥ 0.0884 (closer than third-degree)
related <- subset(related, related$Kinship >= 0.0884)

# Remove individuals who appear in two or more pairs (ID1)
value_counts <- table(related$ID1)
value_to_remove <- as.numeric(names(value_counts[value_counts >= 2]))
baseline <- baseline[!(baseline$eid %in% value_to_remove), ]
related <- related[!(related$ID1 %in% value_to_remove), ]

# Remove individuals who appear in two or more pairs (ID2)
value_counts <- table(related$ID2)
value_to_remove <- as.numeric(names(value_counts[value_counts >= 2]))
baseline <- baseline[!(baseline$eid %in% value_to_remove), ]
related <- related[!(related$ID2 %in% value_to_remove), ]

# Bidirectional pair list (ID1 <-> ID2)
ID1 <- related$ID1
ID2 <- related$ID2
length(intersect(ID1, ID2)) # for reference

tmp <- related
colnames(tmp) <- c("ID2", "ID1", "HetHet", "IBS0", "Kinship")
tmp <- tmp[, c("ID1", "ID2", "HetHet", "IBS0", "Kinship")]
tmp <- rbind(related, tmp)

# Remove remaining participants in ≥2 pairs
value_counts <- table(tmp$ID1)
value_to_remove <- as.numeric(names(value_counts[value_counts >= 2]))
baseline <- baseline[!(baseline$eid %in% value_to_remove), ]
related <- related[!(related$ID1 %in% value_to_remove), ]
rm(tmp)

# Flag who is still present in the baseline dataset
related$ID1_in_baseline <- ifelse(related$ID1 %in% baseline$eid, "Keep", "Remove")
related$ID2_in_baseline <- ifelse(related$ID2 %in% baseline$eid, "Keep", "Remove")

# Drop pairs where both participants were already removed
related <- subset(related, !(ID1_in_baseline == "Remove" & ID2_in_baseline == "Remove"))

# Subset to unresolved related pairs (both in baseline)
related_1 <- subset(related, ID1_in_baseline == "Keep" & ID2_in_baseline == "Keep")
related_1 <- related_1[, c("ID1", "ID2")]

# Add participant metadata
tmp <- baseline[, c("eid", "sex", "enroll_age", "omics")]
related_1 <- merge(related_1, tmp, by.x = "ID1", by.y = "eid")
related_1 <- related_1[, c("ID1", "sex", "enroll_age", "omics", "ID2")]
related_1 <- merge(related_1, tmp, by.x = "ID2", by.y = "eid")
related_1 <- related_1[, c("ID1", "sex.x", "enroll_age.x", "omics.x", "ID2", "sex.y", "enroll_age.y", "omics.y")]

# Prioritize retention based on omics data
table(related_1$omics.x)
table(related_1$omics.y)

rm(tmp)

# Preferentially retain participant with "both" omics
swap_condition <- related_1$omics.y == "both" & related_1$omics.x != "both"
related_1[swap_condition, c("ID1", "ID2")] <- related_1[swap_condition, c("ID2", "ID1")]
related_1[swap_condition, c("omics.x", "omics.y")] <- related_1[swap_condition, c("omics.y", "omics.x")]

# Double-check omics distribution after swap
table(related_1$omics.x)
table(related_1$omics.y)

#~ related_1 has a list of people who are related. we need to figure out which individual from each pair we want to keep
#~ we are going to try and optimize/select by people who have a particular omic set (we can have different DFs for each analysis)

# Remove selected related individuals, prioritizing those with metab data
metab_related <- related_1

# Swap to retain those with "metab" or "both" over "prot" only
swap_condition <- (metab_related$omics.y %in% c("both", "metab")) & (metab_related$omics.x == "prot")
metab_related[swap_condition, c("ID1", "ID2")] <- metab_related[swap_condition, c("ID2", "ID1")]
metab_related[swap_condition, c("omics.x", "omics.y")] <- metab_related[swap_condition, c("omics.y", "omics.x")]

# Final omics distribution (after preference-based swapping)
table(metab_related$omics.x)
table(metab_related$omics.y)

# Drop ID2 individuals to keep one per pair
baseline_metab <- baseline[!(baseline$eid %in% metab_related$ID2), ]

# Remove selected related individuals, prioritizing those with prot data
prot_related <- related_1

# Swap to retain those with "prot" or "both" over "metab" only
swap_condition <- (prot_related$omics.y %in% c("both", "prot")) & (prot_related$omics.x == "metab")
prot_related[swap_condition, c("ID1", "ID2")] <- prot_related[swap_condition, c("ID2", "ID1")]
prot_related[swap_condition, c("omics.x", "omics.y")] <- prot_related[swap_condition, c("omics.y", "omics.x")]

# Final omics distribution (after preference-based swapping)
table(prot_related$omics.x)
table(prot_related$omics.y)

# Drop ID2 individuals to keep one per pair
baseline_prot <- baseline[!(baseline$eid %in% prot_related$ID2), ]

# Clean up environment
rm(related, related_1, metab_related, prot_related, remove, noprev, selected_samples)

###### SHRINK DATAFRAMES ######
baseline_metab <- subset(baseline_metab, baseline_metab$metab == 1)
index <- c(which(colnames(baseline_metab) == "eid"):which(colnames(baseline_metab) == "prot"), which(colnames(baseline_prot) == "omics")) 
baseline_metab <- baseline_metab[, ..index]

baseline_prot <- subset(baseline_prot, baseline_prot$prot == 1) 
index <- c(which(colnames(baseline_prot) == "eid"):which(colnames(baseline_prot) == "metab"),
           which(colnames(baseline_prot) == "double_metab"):which(colnames(baseline_prot) == "omics")) 
baseline_prot <- baseline_prot[, ..index]

###### IMPUTATION AND OUTLIERS ######

### Metabolomics: Imputation + Outlier Removal + Scaling
# Identify column range for metabolomics analytes
start_col <- which(colnames(baseline_metab) == "Clinical_LDL_C")   # 72
end_col   <- which(colnames(baseline_metab) == "Omega_6_pct_PUFA") # 396

# Impute missing values using k-nearest neighbors
metabolomics_imputed <- impute.knn(
  as.matrix(baseline_metab[, start_col:end_col]),
  rowmax = 0.1, colmax = 0.1, k = 10, rng.seed = 1
  )

baseline_metab[, start_col:end_col] <- as.data.frame(metabolomics_imputed$data)

# Remove extreme outliers (±8 SD from median), set them to NA
for (i in start_col:end_col) {
  col <- colnames(baseline_metab)[i]
  values <- baseline_metab[[col]]
  med <- median(values, na.rm = TRUE)
  sd_val <- sd(values, na.rm = TRUE)
  within_bounds <- values > (med - 8 * sd_val) & values < (med + 8 * sd_val)
  baseline_metab[[col]] <- ifelse(within_bounds & !is.na(values), values, NA)
}

# Scale all metabolites (mean = 0, SD = 1)
for (i in start_col:end_col) {
  baseline_metab[, colnames(baseline_metab)[i]] <- scale(as.numeric(as.matrix(baseline_metab[, ..i])))
}

### Proteomics: Scaling Only (no imputation, 20% missingness threshold applied)
# Identify column range for proteins
start_prot <- which(colnames(baseline_prot) == "ERP44")   # 73
end_prot   <- which(colnames(baseline_prot) == "COMMD1")  # 2983

# Scale all proteins
for (i in start_prot:end_prot) {
  baseline_prot[, colnames(baseline_prot)[i]] <- scale(as.numeric(as.matrix(baseline_prot[, ..i])))
}

# Reference checks
which(colnames(baseline_prot) == "ERP44")    # 73
which(colnames(baseline_prot) == "COMMD1")   # 2983

# Clean up
rm(proteomics_imputed, metabolomics_imputed)

###### SAVE OUTPUT FILES ######
# Write metabolomics dataset to compressed TSV
write.table(
  baseline_metab,  # 22,630 participants
  file = gzfile("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.n22630.20240403.tsv.gz"), 
  append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE
  )

# Write proteomics dataset to compressed TSV
write.table(
  baseline_prot,  # 5,023 participants
  file = gzfile("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.noImpute.n5023.20240403.tsv.gz"), 
  append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE
  )

# Final cleanup
rm(b, ID1, ID2, index, swap_condition, value_counts, value_to_remove, save)
