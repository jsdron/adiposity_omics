#########################################
# Script: MR_outputs.R
# Description: Formats and organizes MR results and genetic instruments.
#              Re-adapt the code for the other MR runs.
# Key Outputs:
#   - MR genetic instruments (Supplemental Tables)
#   - Merged MR results (Supplemental Tables)
#########################################

###### LIBRARIES ######
library(data.table)
library(dplyr)

###### SET VARIABLES ######
# Path to MR results
mpan <- "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/results/MR/MET_ADIPOSE/"

mpan_ivs <- c(paste0(mpan, "genetic_instruments_for_exposures_Met-Adipose2024-07-30.csv"))
mpan_results <- c(paste0(mpan, "all_outcomes_mr_results_Met-Adipose2024-07-30.csv"))

###### GENETIC INSTRUMENTS ######
# Read and concatenate IV data
mpan_ivs_data <- do.call(cbind, lapply(mpan_ivs, fread))

# Drop columns with pattern "id.*"
mpan_ivs_data <- mpan_ivs_data[, !grep("^id\\d*$", names(mpan_ivs_data)), with = FALSE]

# Remove duplicated first column
mpan_ivs_data <- mpan_ivs_data[, -1]

# Fix column name ending in .y
colnames(mpan_ivs_data)[1] <- "met-d-Acetate"

# Remove "met-d-" prefix from column names
colnames(mpan_ivs_data) <- gsub("met-d-", "", colnames(mpan_ivs_data))

# Transpose and merge each row into comma-separated strings
transposed <- as.data.frame(t(mpan_ivs_data))
transposed <- apply(transposed, 1, function(row) {
  row <- row[!is.na(row) & row != ""]
  paste(row, collapse = ", ")
})
transposed <- data.frame(transposed, stringsAsFactors = FALSE)

# Split into MR and MR-RAPS instrument sets
mr_iv <- transposed[-grep("MR.RAPS", rownames(transposed)), , drop = FALSE]
mr_rap_iv <- transposed[grep("MR.RAPS", rownames(transposed)), , drop = FALSE]

# Remove duplicates
mr_iv <- mr_iv[!duplicated(mr_iv), , drop = FALSE]
mr_rap_iv <- mr_rap_iv[!duplicated(mr_rap_iv), , drop = FALSE]

# Save instrument tables
write.table(mr_iv, "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/results/MR/Instruments/all_ivs_7.30_Metab.table.tsv",
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

write.table(mr_rap_iv, "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/results/MR/Instruments/all_raps_ivs_7.30_Metab.table.tsv",
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")


###### MR RESULTS ######
# Read and concatenate MR results
mpan_results_data <- do.call(rbind, lapply(mpan_results, fread))

# Save merged MR results
write.table(mpan_results_data, "/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/results/MR/Instruments/all_merged_MR_results_7.30_Metab.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
