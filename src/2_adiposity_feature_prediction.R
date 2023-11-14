### 0 - LIBRARIES ####

  library(data.table)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)

### 1 - LOAD AND MERGE DATAFRAMES ####
  
  omic <- c("metabol", "proteom")
  i = 1
  
  df <-fread(paste0("../input/adiposity_",omic[i],"_df.csv"))

### 2 - SET VARIABLE LISTS ####

  met_features <- colnames(df)[which(colnames(df) == "Clinical_LDL_C"):which(colnames(df) == "Omega_6_pct_PUFA")]
  # prot_features
  
  j = 1

### 3 - FEATURE ASSOCIATION ####

  adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height
  
  k = 1

  model <- lm(adi_traits[k] ~)


# Loop through lagou and fit GLM models
for (j in omic_features) {
  formula <-
    as.formula(paste0("has_obesity ~", var, "+ enroll_age + mri + time between + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
  model <-
    glm(formula, data = file, family = binomial) # Change family and response_variable as needed
  output <-
    tidy(
      model,
      conf.int = FALSE,
      conf.level = 0.95,
      exponentiate = TRUE
    )
  result_row <- output[which(output$term == var), ]
  
  # Append the row to result_df
  glm_Obs <- rbind(glm_Obs, result_row)
  
}




ggplot(joined_df, aes(x = OR_Obs, y = OR_CAD, color = sig)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "black",
    linetype = "dashed"
  ) +
  # scale_shape_manual(values=c(1,19)) +
  geom_errorbar(aes(ymin = lower_ci_Obs, ymax = upper_ci_Obs),
                width = 0.1,
                color = "#F2F2F2") +  # Vertical error bars
  geom_errorbarh(aes(xmin = lower_ci_CAD, xmax = upper_ci_CAD),
                 height = 0.05,
                 color = "#F2F2F2") +  # Horizontal error bars
  
  # Make the axes look better
  #   scale_y_continuous(limits = c(1, 3), breaks = seq(0, 3, by =1)) +
  #   scale_x_continuous(limits = c(1, 7), breaks = seq(0, 7, by = 1)) +
  
  geom_point(aes(
    x = OR_Obs,
    y = OR_CAD,
    color = sig,
    size = 3
  )) +
  scale_color_manual(values = c("#FF5000", "#006DB6", "#F2F2F2")) +
  
  # Plot the ones I want on top!
  geom_point(data = filter(joined_df, sig != "NS"),
             aes(
               x = OR_Obs,
               y = OR_CAD,
               color = sig,
               size = 3
             )) +
  
  # Hide legends I don't want
  guides(size = FALSE, color = FALSE) +
  
  # Add labels
  # geom_text(aes(x = OR_Obs, y = OR_CAD, label = label), vjust = -0.5, hjust = -0.5, color = "black") +
  
  geom_text_repel(
    data = filter(joined_df, any_sig == "Significant"),
    aes(label = label),
    vjust = -0.5,
    hjust = -0.5,
    color = "black",
    box.padding = unit(0.2, "lines")
  ) +
  
  coord_cartesian(ylim = c(1, 3), xlim = c(1, 7)) +  # Set axis limits
  
  # Add labels for axes
  labs(
    x = "Obesity\nOdds Ratio (95% CI)",
    y = "Coronary Artery Disease\nOdds Ratio (95% CI)",
    color = "Significance",
    shape = "Significance"
  ) +
  
  #theme(legend.position = "bottom", legend.justification = "right") +
  theme_cowplot()