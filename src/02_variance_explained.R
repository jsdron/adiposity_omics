#########################################
# Script: variance_explained.R
# Description: This script reads in proteomics data,
#              prepares the input matrix, and computes the number
#              of principal components required to explain 90% and 99% of the variance.
# Key Output: Printed number of components explaining 90% and 99% variance --> can use to determine correction factor for Bonferroni significance
#########################################

###### LIBRARIES ######
library(data.table)
library(tidyr)
library(dplyr)

###### FUNCTIONS ######
# Function to calculate and report the number of PCA components that explain 90% and 99% of the variance
pca_variance_explained <- function(matrix) {
  # Perform PCA
  pca_result <- prcomp(matrix, scale. = TRUE)
  
  # Calculate the cumulative variance explained by the components
  cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
  
  # Find the number of components explaining 99% variance (or 90%)
  num_components_99 <- which(cumulative_variance >= 0.99)[1]
  num_components_90 <- which(cumulative_variance >= 0.90)[1]
  
  # Print the number of components explaining 99% and 90% variance
  cat("Number of components explaining 99% variance:", num_components_99, "\n")
  cat("Number of components explaining 90% variance:", num_components_90, "\n")
}

###### PREPARE DATAFRAME(S) ######
# Load baseline proteomics data
baseline_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.noImpute.n5023.20240403.tsv.gz') 

# Standardize column names
colnames(baseline_prot) <- gsub("-", "_", colnames(baseline_prot))

# Select the first and last analyte columns
analyte <- colnames(baseline_prot)[which(colnames(baseline_prot) == "ERP44"):which(colnames(baseline_prot) == "COMMD1")]

# Convert to matrix and remove rows with missing values
matrix <- as.matrix(baseline_prot[, ..analyte])
matrix_clean <- matrix[complete.cases(matrix), ]

###### VARIANCE EXPLAINED VIA PCA ######
pca_variance_explained(matrix_clean)