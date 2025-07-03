#########################################
# Script: table1.R
# Description: This script combines metabolomic and proteomic baseline data,
#              adds labels for diabetes and waist-to-hip ratio (WHR), and generates a summary Table 1
#              stratified by omics group (metabolomics vs. proteomics).
# Key Output: Table 1 (stratified summary statistics)
#########################################

###### LIBRARIES ######
library(data.table)
library(tableone)

###### PREPARE DATAFRAME(S) ######
# Load baseline data for the metabolomic and proteomic subsets
baseline_metab <- fread("/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/data/baseline_adi_metab.n22630.20240403.tsv.gz")
baseline_prot <- fread("/Users/michael.tian/Desktop/Natarajan_Lab/adiposity/adiposity_omics/data/baseline_adi_prot.noImpute.n5023.20240730.tsv.gz")

# Add sex label
baseline_metab$sex_label <- ifelse(baseline_metab$sex == 1,"Male",
                                   ifelse(baseline_metab$sex == 0, "Female", NA))
baseline_prot$sex_label <- ifelse(baseline_prot$sex == 1,"Male",
                                  ifelse(baseline_prot$sex == 0, "Female", NA))

# Add omic group label
baseline_metab$trt = "met"
baseline_prot$trt = "prot"

# Combine datasets
Overall_data = rbind(baseline_metab, baseline_prot, fill = TRUE)

# Adding diabetes information
hba1c <- fread("/Users/michael.tian/Desktop/Natarajan_Lab/Adiposity/Adiposity_Omics/data/UKBB_TG_BackgroundVariables_22020.txt",
    select = c('eid', 'HgbA1c'))

# Combine datasets
Overall_data <- merge(Overall_data, hba1c, by = "eid")


# Define diabetes and prediabetes labels (38.8 is the cutoff for prediabetes and 47.5 is the cutoff for diabetes)
Overall_data$Diabetes <- ifelse(Overall_data$HgbA1c < 47.5, "Normal", "Diabetes")
Overall_data$Prediabetes <- ifelse(Overall_data$HgbA1c < 38.8, "Normal",
    ifelse(Overall_data$HgbA1c < 47.5, "Pre-Diabetes", "Normal"))

# Define WHR status
Overall_data$High_WHR_Ratio <- ifelse(Overall_data$WHRInstance2 > 0.85, "High", "Normal")
Overall_data$High_WHR_Ratio[Overall_data$WHRInstance2 <= 0.9 & Overall_data$sex_label == "Male"] <- "Normal"

###### GENERATE TABLE 1 ######
# Set factor levels for ordered displays
Overall_data$sex_label <- factor(Overall_data$sex_label, levels = c("Male", "Female"))
Overall_data$Diabetes <- factor(Overall_data$Diabetes, levels = c("Normal", "Diabetes"))
Overall_data$Prediabetes <- factor(Overall_data$Prediabetes, levels = c("Normal", "Pre-Diabetes"))
Overall_data$High_WHR_Ratio <- factor(Overall_data$High_WHR_Ratio, levels = c("Normal", "High"))

# Define variables to include in Table 1
myVars <- c(
    "enroll_age",
    "sex_label",
    'currentsmoker',
    "fasting_status",
    "BMI",
    "BMI_class",
    "time_between",
    "Diabetes",
    "Prediabetes",
    "WHRInstance2",
    "High_WHR_Ratio",
    "asatadjbmi",
    "gfatadjbmi",
    "vatadjbmi"
    )

# Define categorical variables
catVars <- c(
    "sex_label",
    'currentsmoker',
    "fasting_status",
    "BMI_class",
    "Diabetes",
    "Prediabetes",
    "High_WHR_Ratio"
    )

# Generate Table 1
tab1 <- CreateTableOne(vars = myVars, data = Overall_data, strata = "trt", factorVars = catVars, addOverall = FALSE)

# Final formatted Table 1 output
table1 <- print(tab1,
        quote = FALSE,
        noSpaces = TRUE,
        printToggle = FALSE #nonnormal = nonnorm,
        )

# Display Table 1
table1
