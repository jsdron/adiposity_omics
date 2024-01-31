#### 0 - LIBRARIES ####

# only need to do these once. comment them out and save after they have been successfully run
install.packages("OneSampleMR")
update.packages()

# run these every time
library("OneSampleMR")
library(data.table)
library(tidyr)
library(dplyr)

#### 1 - LOAD AND MERGE DATAFRAMES ####

### Baseline covariates
baseline <- fread("/Volumes/medpop_esp2/aniruddh/Lpa/UKBB_TG_BackgroundVariables_22020.txt",
                  select=c('eid','enroll_age','sex','ethnicity','eversmoker','eversmoked','currentsmoker', 
                           'DBP','SBP','fasting1','anylipidmed0','cholmed0','statin0'))

tmp <- fread("/Volumes/medpop_esp2/aschuerm/resources/ukb_df_base_full_20230905.tsv.gz",  
             select=c("id", "BMI", "race", "antihtnbase", "BP_Meds"))

baseline <- merge(baseline, tmp, by.x="eid", by.y="id")
rm(tmp)

baseline$fasting_status <- ifelse(baseline$fasting1 >=8, 1, 0) # considered fasted (1) if fasting time is greater than or equal to 8 hours
baseline$BMI_class <- ifelse(baseline$BMI < 18.5, "underweight",
                             ifelse(baseline$BMI >= 18.5 & baseline$BMI < 25, "normal_weight",
                                    ifelse(baseline$BMI >= 25 & baseline$BMI < 30, "overweight", "obesity")))

colnames(baseline) <- c("eid","enroll_age","sex","ethnicity", "eversmoker","eversmoked","currentsmoker",
                        "DBP", "SBP", "fasting_time",'anylipidmed','chol_med','statin','BMI','race',
                        'HTN_med','BP_med', 'fasting_status','BMI_class')

baseline <- baseline[,c("eid","enroll_age","sex","ethnicity",'race',"eversmoker","eversmoked","currentsmoker",
                        'BMI','BMI_class',"DBP", "SBP", "fasting_time", 'fasting_status','anylipidmed','chol_med','statin',
                        'HTN_med','BP_med')]

### MRI-derived adiposity measurements
adi <- fread('/Volumes/medpop_esp2/sagrawal/shared_fat/data/fatPhenoCovar_wFMR.csv') # 40032
adi <- adi[,-c("sex")]

# exclude people w/o adiposity measurements
baseline <- subset(baseline, !(baseline$eid %in% adi$eid))
rm(adi)

#### 2 - REMOVE WITHDRAWN ####
withdrawn <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/withdrawn_samples/w7089_2023-08-21.csv')
colnames(withdrawn)[1] <- "eid"
withdrawn_vector <- c(withdrawn$eid)
baseline <- baseline[!(baseline$eid %in% withdrawn_vector),] # 40021, 11 withdrew
rm(withdrawn, withdrawn_vector)

#### 3 - KEEP EUR ONLY ####
baseline <- subset(baseline, baseline$race=="white") # 38701; 1320 removed

#### 4 - KEEP PARITICIPANTS WITH PROT OR METAB MEASUREMENTS ####

### Metabolomics data 
qc_met <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb49184/ukbnmr/ukb49184_metabolome_qced.tsv.gz') # QC'd metabolites from Satoshi, using ukbnmr() package to remove technical artifacts
double_measures <- qc_met$eid[duplicated(qc_met$eid)]
qc_met <- subset(qc_met, qc_met$visit_index==0) # take the baseline measurement 
index <- c(which(colnames(qc_met)=="eid"), 
           which(colnames(qc_met)=="Clinical_LDL_C"):which(colnames(qc_met)=="Omega_6_pct_PUFA")) # create a dataframe (i.e., b) that only contains id's and metabolomics data
b <- qc_met[, ..index]
b <- b %>% select(where(~mean(is.na(.)) <= 0.1)) # remove metabolites with more than 10% missingness
b <- b[rowMeans(is.na(b[,-1])) <= 0.1,] # exclude individuals with more than 10% missing metabolites ( ppl removed)
qc_met <- qc_met[qc_met$eid %in% b$eid,]

baseline$metab <- ifelse(baseline$eid %in% qc_met$eid, 1, 0) # 1 if there is metab data available for the participant
baseline <- merge(baseline, qc_met, by="eid", all.x=TRUE)
baseline$double_metab <- ifelse(baseline$eid %in% double_measures, 1, 0) # if participant has 2 metabolomic measurements, indicate with a 1
rm(qc_met, double_measures, index, b)

### Proteomics data ([1] load data, [2] load linker file to assign protein names, [3] select only those with baseline protein measurements, [4] make dataframe wide, [5] rename columns)
prot <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb672544/olink_data.txt.gz')
names_prot <- fread('/Volumes/medpop_esp2/aschuerm/resources/ukb/coding_ukb_olink.tsv')
prot <- prot[prot$ins_index==0, -("ins_index")]
prot <- as.data.frame(pivot_wider(prot, names_from = protein_id, values_from = result))
names_prot$coding <- as.character(names_prot$coding)
setnames(prot, old = names_prot$coding, new = names_prot$short_name, skip_absent = TRUE)

index <- c(which(colnames(prot)=="eid"),
           which(colnames(prot)=="CLIP2"):which(colnames(prot)=="NPM1")) # create a dataframe (i.e., b) that only contains id's and proteomics data
b <- prot[, index]
c <- b %>% select(where(~mean(is.na(.)) <= 0.1)) # remove proteins with more than 10% missingness
setdiff(names(b), names(c)) # proteins removed: PCOLCE, TACSTD2, CTSS, NPM1
c <- c[rowMeans(is.na(b[,-1])) <= 0.1,] # exclude individuals with more than 10% missing proteins (N=4051)
c <- c[, !names(c) %in% c('PCOLCE', 'TACSTD2', 'CTSS', 'NPM1')]

baseline$prot <- ifelse(baseline$eid %in% c$eid, 1, 0) # 1 if there is proteomic data available for the participant
baseline <- merge(baseline, c, by="eid", all.x=TRUE)

# individuals with either prot or metab data
baseline <- subset(baseline, baseline$prot==1 | baseline$metab==1) # 
baseline$omics <- ifelse(baseline$metab==1 & baseline$prot==1, "both", 
                         ifelse(baseline$metab==1 & baseline$prot==0, "metab", "prot"))

rm(prot, names_prot, index, b, c)

#### 5 - REMOVE RELATIVES ####
related <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb7089_rel_s488264.dat') # Mesbah passed me this path

# Following what I used in MGBB for GLGC, I am removing related individuals with a kinship >=0.0884 (anything closer than 3-degree relative)
related <- subset(related, related$Kinship>=0.0884) # Each line in the related df = a pair or related individuals

value_counts <- table(related$ID1) # Count occurrences of each unique value
value_to_remove <- as.numeric(names(value_counts[value_counts >= 2])) # defs exclude these individuals b/c they are related to 2 or more individuals
baseline <- baseline[!(baseline$eid %in% value_to_remove),] # exclude the ppl related to 2 or more individuals 
related <- related[!(related$ID1 %in% value_to_remove),]

value_counts <- table(related$ID2) 
value_to_remove <- as.numeric(names(value_counts[value_counts >= 2])) # defs exclude these individuals b/c they are related to 2 or more individuals
baseline <- baseline[!(baseline$eid %in% value_to_remove),] # exclude the ppl related to 2 or more individuals (N=12591)
related <- related[!(related$ID2 %in% value_to_remove),]

ID1 <- related$ID1
ID2 <- related$ID2
length(intersect(ID1, ID2)) # 537

tmp <- related
colnames(tmp) <- c("ID2", "ID1","HetHet","IBS0","Kinship")
tmp <- tmp[,c("ID1", "ID2","HetHet","IBS0","Kinship")]
tmp <- rbind(related, tmp)

value_counts <- table(tmp$ID1) # Count occurrences of each unique value
value_to_remove <- as.numeric(names(value_counts[value_counts >= 2])) # defs exclude these individuals b/c they are related to 2 or more individuals
baseline <- baseline[!(baseline$eid %in% value_to_remove),] # exclude the ppl related to 2 or more individuals (N=12575)
related <- related[!(related$ID1 %in% value_to_remove),]
rm(tmp)

# There are 30595 related pairs... to figure out which single person of each pair to keep:
# A) exclude those who are NOT in the current baseline dataframe 
related$ID1_in_baseline <- ifelse(related$ID1 %in% baseline$eid, "Keep", "Remove")
related$ID2_in_baseline <- ifelse(related$ID2 %in% baseline$eid, "Keep", "Remove")
# If both columns are "Remove", then remove the line
related <- subset(related, !(ID1_in_baseline == "Remove" & ID2_in_baseline == "Remove"))
# For all the pairs where 1 is Keep and 1 is Remove, I don't have to remove anything from the main dataframe, because it is already restricted to 1 participant per pair
# For the pairs where both are listed as Keep, I need to choose 1 participant to keep
related_1 <- subset(related, ID1_in_baseline == "Keep" & ID2_in_baseline == "Keep")
related_1 <- related_1[,c("ID1", "ID2")]

tmp <- baseline[,c("eid","sex","enroll_age","omics")]
related_1 <- merge(related_1, tmp, by.x="ID1", by.y="eid")
related_1 <- related_1[,c("ID1","sex","enroll_age","omics","ID2")]
related_1 <- merge(related_1, tmp, by.x="ID2", by.y="eid")
related_1 <- related_1[,c("ID1","sex.x","enroll_age.x","omics.x","ID2","sex.y","enroll_age.y","omics.y")]

table(related_1$omics.x)
table(related_1$omics.y)

rm(tmp)

# Create a condition for swapping based on "both" in omics.y
swap_condition <- related_1$omics.y == "both" & related_1$omics.x != "both"

# Swap values for ID1 and ID2, and omics.x and omics.y where the condition is met
related_1[swap_condition, c("ID1", "ID2")] <- related_1[swap_condition, c("ID2", "ID1")]
related_1[swap_condition, c("omics.x", "omics.y")] <- related_1[swap_condition, c("omics.y", "omics.x")]

table(related_1$omics.x)
table(related_1$omics.y)

## Metab
metab_related <- related_1 
swap_condition <- (metab_related$omics.y == "both" | metab_related$omics.y == "metab") & (metab_related$omics.x =="prot")

# Swap values for ID1 and ID2, and omics.x and omics.y where the condition is met
metab_related[swap_condition, c("ID1", "ID2")] <- metab_related[swap_condition, c("ID2", "ID1")]
metab_related[swap_condition, c("omics.x", "omics.y")] <- metab_related[swap_condition, c("omics.y", "omics.x")]

table(metab_related$omics.x)
table(metab_related$omics.y)

# ID1 now has most with both prot or both data
baseline_metab <- baseline[!(baseline$eid %in% metab_related$ID2),] # exclude ID2

## Prot
prot_related <- related_1 
swap_condition <- (prot_related$omics.y == "both" | prot_related$omics.y == "prot") & (prot_related$omics.x =="metab")

# Swap values for ID1 and ID2, and omics.x and omics.y where the condition is met
prot_related[swap_condition, c("ID1", "ID2")] <- prot_related[swap_condition, c("ID2", "ID1")]
prot_related[swap_condition, c("omics.x", "omics.y")] <- prot_related[swap_condition, c("omics.y", "omics.x")]

table(prot_related$omics.x)
table(prot_related$omics.y)

# ID1 now has most with both prot or both data
baseline_prot <- baseline[!(baseline$eid %in% prot_related$ID2),] # exclude ID2

rm(ID1, ID2, swap_condition, value_counts, value_to_remove)

#### 6 - FINALIZE FRAMES ####
# Metabolomics
rm(baseline)
baseline_metab <- subset(baseline_metab, baseline_metab$metab==1)

index <- c(which(colnames(baseline_metab)=="eid"):which(colnames(baseline_metab)=="visit_index"), 
           which(colnames(baseline_metab)=="Clinical_LDL_C"):which(colnames(baseline_metab)=="Omega_6_pct_PUFA")) # create a dataframe (i.e., b) that only contains id's and metabolomics data
baseline_metab <- baseline_metab[, ..index]

write.table(baseline_metab, # 97101
            file = gzfile("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_1MR.n97101.20240116.tsv.gz"), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

# Proteomics
baseline_prot <- subset(baseline_prot, baseline_prot$prot==1)
index <- c(which(colnames(baseline_prot)=="eid"):which(colnames(baseline_prot)=="BP_med"), 
           which(colnames(baseline_prot)=="prot"):which(colnames(baseline_prot)=="SCARB2")) # create a dataframe (i.e., b) that only contains id's and metabolomics data
baseline_prot <- baseline_prot[, ..index]

write.table(baseline_prot, # 4684
            file = gzfile("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/prot_1MR.n40015.20240116.tsv.gz"), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


#### 7 - GWAS FILES ####

### Functions ###
# Code from Saaket to run inverse-normalize rank transformation of the data
rINTtransform <- function(vec){
  y<-qnorm((rank(vec,na.last='keep')-0.5)/sum(!is.na(vec)))
  return(y)
}

## takes the residuals from the model and then inverse-rank normalizes them
transformAnalyte <- function(analyte_name, data){
  formula <- as.formula(paste(analyte_name, "~ enroll_age + sex + enroll_age2 +  
                               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  
  mdl <- glm(formula, data = data, na.action = na.omit)
  x <- resid(mdl)
  
  return(rINTtransform(x)) 
}

## PCs
pcs <- fread("/Volumes/medpop_esp2/skoyama/passing/ukb_kgpprojection/v02/out/ukb.kgp_projected.tsv.gz")

### Metabolomics (Nightengale) ###
    # https://www.nature.com/articles/s41586-023-05844-9#Sec10 
    metab <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_1MR.n97101.20240116.tsv.gz")
    which(colnames(metab)=="Clinical_LDL_C") # The first metabolite: column 22
    which(colnames(metab)=="Omega_6_pct_PUFA") # The last metabolite: column 346
    
    ## Remove extreme outliers +/- 8 SD from the median by setting them to NA.
    for (i in 22:346) {
      col <- colnames(metab)[i]
      metab[[col]] <- ifelse(metab[[col]] > (median(metab[[col]], na.rm=TRUE) - 8*sd(metab[[col]], na.rm=TRUE)) & 
                               metab[[col]] < (median(metab[[col]], na.rm=TRUE) + 8*sd(metab[[col]], na.rm=TRUE)) & 
                                !is.na(metab[[col]]), metab[[col]], NA)
    }
    
    ## Natural-log transform each metabolite
    which(colnames(metab)=="Clinical_LDL_C") # The first metabolite: column 22
    which(colnames(metab)=="Remnant_C") # Before the percentages start: column 189
    which(colnames(metab)=="Omega_6_pct_PUFA") # The last metabolite: column 346
    
    norm_file <- metab %>% as_tibble() %>% mutate(across(c(22:189), ~log(. + 1e-10))) # Adding a small constant to avoid issues with zero or negative values. This code adds a small constant (1e-10) to each value before taking the logarithm. Adjust the constant based on your specific needs.
    norm_file <- merge(norm_file, pcs, by="eid") # Add PCs
    norm_file$enroll_age2 <- norm_file$enroll_age^2
    
    transformed_results <- norm_file
    
    ## Generate residuals
    for (column_name in names(transformed_results)[22:189]) {
      # checking for NAs
      original_column <- norm_file[[column_name]]
      non_na_indices <- which(!is.na(original_column))
      
      transformed_values <- transformAnalyte(column_name, data = transformed_results)
      # Align transformed values with original column based on non-missing values
      transformed_results[[column_name]][non_na_indices] <- transformed_values
      
    }
    
    ## Create file for PHENOTYPES that will be run through GWAS
    FID <- 0
    gwas_df <- data.frame(FID, transformed_results[,c(1,22:189)])
    gwas_df$FID <- 0
    colnames(gwas_df)[2] <- "IID"
    gwas_df$FID <- gwas_df$IID
    
    # Save as a new file that will be used specifically for GWAS
    write.table(gwas_df, # 96613
                file = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_gwas.PHENO.n96613.20240123.txt", 
                append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
    
    # Create list of the phenotypes that need to be looped through for GWAS
    list <- colnames(gwas_df)
    list <- list[3:length(list)]
    write.table(list, # 96613
                file = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_gwas.list.txt", 
                append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
  
    ## Create file for COVARIATE that will be run through GWAS
    FID <- 0
    gwas_df <- data.frame(FID, 
                          transformed_results[,c("eid","sex","enroll_age","enroll_age2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","fasting_time","currentsmoker")]) 
    gwas_df$FID <- 0
    colnames(gwas_df)[2] <- "IID"
    gwas_df$FID <- gwas_df$IID
    
    # Add additional covariates
    baseline <- fread("/Volumes/medpop_esp2/aniruddh/Lpa/UKBB_TG_BackgroundVariables_22020.txt",
                      select=c('eid','cholmed0','fishoil0')) #statin0
        
    colnames(baseline) <- c("eid","cholmed","fishoil")
    gwas_df <- merge(gwas_df, baseline, by.x="FID", by.y="eid")
 
    write.table(gwas_df, # 96613
                file = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_gwas.COVAR.n96613.20240123.tsv", 
                append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
    
    

### Proteomics (Olink) ###
# transform the variables!!! use the same code as above, but modified for the protein file
# then write the transformed dataframe to a NEW file. do not overwrite the old one!!
