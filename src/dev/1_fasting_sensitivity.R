### 0 - LIBRARIES ####

  library(data.table)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(broom)

### 1 - LOOP THROUGH LINEAR REGRESSION MODELS ####

  omic <- c("metabol", "proteom")
  result_row <- data.frame()
  
  for (i in 1:2) {
    df <- fread(paste0("../input/adiposity_", omic[i], "_df.csv"))
    
    if (omic[i] == "metabol") {
      analyte <- colnames(df)[which(colnames(df) == "Clinical_LDL_C"):which(colnames(df) == "Omega_6_pct_PUFA")] # metab
    } else {
      colnames(df) <- gsub("-", "_", colnames(df))
      analyte <- colnames(df)[which(colnames(df) == "CLIP2"):which(colnames(df) == "SCARB2")] # prot
    }
    
    for (j in 1:length(analyte)) {
      adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3")
      
      for (k in 1:3) {
        formulas <- list(
          base_formula = paste0(adi_traits[k], " ~ ", analyte[j], "+ sex + age + time_between + mriNum + race + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
          fast_stat_formula = paste0(adi_traits[k], " ~ ", analyte[j], "+ sex + age + time_between + mriNum + race + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + fasting_status"),
          fast_time_formula = paste0(adi_traits[k], " ~ ", analyte[j], "+ sex + age + time_between + mriNum + race + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + fasting_time")
        )
        
        for (formula_name in names(formulas)) {
          model <- lm(as.formula(formulas[[formula_name]]), data = df)
          output <- tidy(model, conf.int = TRUE, conf.level = 0.95)
          output$trait <- adi_traits[k]
          output$model <- formula_name
          tmp <- output[which(output$term == analyte[j]), ]
          result_row <- rbind(result_row, tmp)
        }
      }
    }
  }

  # write.csv(result_row, "../data/fasting_effects_allAdipose_allOmics.csv", row.names=F)

  
### 2 - VISUALIZE AND EVALUATE RESULTS ####  
  