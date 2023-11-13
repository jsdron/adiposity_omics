### 0 - LIBRARIES ####

  library(data.table)
  library(tidyr)
  library(dplyr)

### 1 - LOAD AND MERGE DATAFRAMES ####

  ### Adipose data 
    adi <- fread('/medpop/esp2/sagrawal/shared_fat/data/fatPhenoCovar_wFMR.csv')

  ### Baseline covariates
    baseline <- fread("/medpop/esp2/aschuerm/resources/ukb_df_base_full_20230905.tsv.gz")
    baseline <- baseline[,c("id", "age", "ever_smoked", "BMI", "SBP", "DBP", 
                            "hdl", "tchol", "tdi", "cholmed", "race", 
                            "antihtnbase", "dm2_prev", "mort_fu", "mort_inc", 
                            "cad_prev", "cad_inc", "cad_fu")]
    adi <- merge(adi, baseline, by.x="eid", by.y="id")
    
  ### Metabolomics data (generating a variable to go along with metabolomic data, indicating who has the data available [1] and who doesn't [0])
    qc_met <- fread('/medpop/esp2/projects/UK_Biobank/baskets/ukb49184/ukbnmr/ukb49184_metabolome_qced.tsv.gz')
    qc_met <- subset(qc_met, qc_met$visit_index==0) # take the initial measurement 
    adi_met <- merge(adi, qc_met, by="eid")

  ### Proteomics data ([1] load data, [2] load linker file to assign protein names, [3] select only those with baseline protein measurements, [4] make dataframe wide, [5] rename columns)
    prot <- fread('/medpop/esp2/projects/UK_Biobank/baskets/ukb672544/olink_data.txt.gz')
    names_prot <- fread('/medpop/esp2/aschuerm/resources/ukb/coding_ukb_olink.tsv')
    prot <- prot[prot$ins_index==0, -("ins_index")]
    prot <- as.data.frame(pivot_wider(prot, names_from = protein_id, values_from = result))
    names_prot$coding <- as.character(names_prot$coding)
    setnames(prot, old = names_prot$coding, new = names_prot$short_name, skip_absent = TRUE)
    adi_prot <- merge(adi, prot, by="eid")
    
    rm(prot, names_prot, qc_met, adi, baseline)

### 2 - EXCLUSIONS ####

  ### Metabolomics dataframe
    dim(adi_met)  # N=9,654
    adi_met <- adi_met[!is.na(adi_met$race),]
    dim(adi_met)  # N=9,634
    adi_met <- adi_met[!is.na(adi_met$PC1),]
    
    dim(adi_met)  # N=9,602
    
    index <- c(which(colnames(adi_met)=="eid"), which(colnames(adi_met)=="Clinical_LDL_C"):which(colnames(adi_met)=="Omega_6_pct_PUFA")) # create a dataframe (i.e., b) that only contains id's and metabolomics data
    b <- adi_met[, ..index]
    b <- b %>% select(where(~mean(is.na(.)) <= 0.1)) # remove metabolites with more than 10% missingness (0 metabolites removed)
    b <- b[rowMeans(is.na(b[,-1])) <= 0.1,] # exclude individuals with more than 10% missing metabolites (N=9,382)
    rm(index)
    
    adi_met <- adi_met[adi_met$eid %in% b$eid,]
    dim(adi_met)  # N=9,382
    
    
  ### Proteomics dataframe
    dim(adi_prot)  # N=5,228
    adi_prot <- adi_prot[!is.na(adi_prot$race),]
    dim(adi_prot)  # N=5,218
    adi_prot <- adi_prot[!is.na(adi_prot$PC1),]
    dim(adi_prot)  # N=5,167
    
    index <- c(which(colnames(adi_prot)=="eid"), which(colnames(adi_prot)=="CLIP2"):which(colnames(adi_prot)=="NPM1")) # create a dataframe (i.e., b) that only contains id's and metabolomics data
    b <- adi_prot[, ..index]
    b <- b %>% select(where(~mean(is.na(.)) <= 0.1)) # remove proteins with more than 10% missingness (4 proteins removed: PCOLCE, TACSTD2, CTSS, NPM1)
    b <- b[rowMeans(is.na(b[,-1])) <= 0.1,] # exclude individuals with more than 10% missing metabolites (N=4,844)
    rm(index)
    
    adi_prot <- adi_prot[adi_prot$eid %in% b$eid,]
    adi_prot <- adi_prot[,-c("PCOLCE", "TACSTD2", "CTSS", "NPM1")]
    dim(adi_prot)  # N=4,844
    

### 3 - IMPUTATION OF METABOLITES AND PROTEINS ####
    
  ### Metabolomics dataframe
    # This is KNN-based imputation (K=10)
    library(impute)
    metabolomics_imputed <- impute.knn(as.matrix(adi_met[,which(colnames(adi_met)=="Clinical_LDL_C"):which(colnames(adi_met)=="Omega_6_pct_PUFA")]), rowmax = 0.1, colmax = 0.1, k = 10, rng.seed=1)
    adi_met[,which(colnames(adi_met)=="Clinical_LDL_C"):which(colnames(adi_met)=="Omega_6_pct_PUFA")] <- as.data.frame(metabolomics_imputed$data)
    # this just scales all metabolites (mean=0 / SD=1) so we can use clean variables in our analyses
    for (i in (which(colnames(adi_met)=="Clinical_LDL_C"):which(colnames(adi_met)=="Omega_6_pct_PUFA"))) {
     adi_met[,colnames(adi_met)[i]] <- scale(as.numeric(as.matrix(adi_met[,..i])))
    }

  ### Proteomics dataframe
    # This is KNN-based imputation (K=10)
    proteomics_imputed <- impute.knn(as.matrix(adi_prot[,which(colnames(adi_prot)=="CLIP2"):which(colnames(adi_prot)=="SCARB2")]), rowmax = 0.1, colmax = 0.1, k = 10, rng.seed=1)
    adi_prot[,which(colnames(adi_prot)=="CLIP2"):which(colnames(adi_prot)=="SCARB2")] <- as.data.frame(proteomics_imputed$data)
    # this just scales all proteins (mean=0 / SD=1) so we can use clean variables in our analyses
    for (i in (which(colnames(adi_prot)=="CLIP2"):which(colnames(adi_prot)=="SCARB2"))) {
     adi_prot[,colnames(adi_prot)[i]] <- scale(as.numeric(as.matrix(adi_prot[,..i])))
    }

    
  ### Write CSV's
    write.csv(adi_met, "/medpop/esp2/aschuerm/ukb_adiposity_omics/input/adiposity_metabol_df.csv", row.names=F)
    write.csv(adi_prot, "/medpop/esp2/aschuerm/ukb_adiposity_omics/input/adiposity_proteom_df.csv", row.names=F)

    