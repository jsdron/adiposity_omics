#### 0 - LIBRARIES ####

library(data.table)
library(tidyr)
library(dplyr)

#### 1 - LOAD AND MERGE DATAFRAMES ####
# Get a list of all the sample IDs that were removed from the genotype/sample QC from scripts 01a-c
keep_samples <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_gwas.PHENO.n96613.20240123.tsv', 
                      select=c("IID"))

#### 2 - REMOVE BAD IDS ####
chr <- c(1:22, "X")

for (i in chr){
  tmp <- fread(paste0('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/tmp/ukb.tmp.chr',i,'.mindrem.id'))
  
  keep_samples <- subset(keep_samples, !(keep_samples$IID %in% tmp$V2))
  
}

#### 3 - FORMAT AND SAVE ####
keep_samples$ID_2 <- keep_samples$IID
colnames(keep_samples)[1] <- "ID_1"

write.table(keep_samples, # 96613
            file = "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/metab_gwas.keep_samples.n82680.20240123.txt", 
            append = FALSE, sep = "\t", dec = ".", col.names = FALSE, row.names = FALSE, quote = FALSE)
