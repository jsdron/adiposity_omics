#### 0 - LIBRARIES ####

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tableone)
library(stringr)
library(broom)
library(ggpubr)

#### 1 - LOAD DATAFRAMES ####

baseline_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.n4684.20231201.tsv.gz') # N = 4684
colnames(baseline_prot) <- gsub("-", "_", colnames(baseline_prot))

#### 2 - SET VARIABLE LISTS ####

analyte <- colnames(baseline_prot)[which(colnames(baseline_prot) == "CLIP2"):which(colnames(baseline_prot) == "SCARB2")]

adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

threshold = 41+1459 # the number of proteins being considered and 41 from metabs


#### 3 - PCs or no PCs ####
## PCs
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum + time_between + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    model <- lm(lm_formula, data = baseline_prot) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}

output_prot$adi_label <- ifelse(output_prot$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

## No PCs
output_prot2 <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_prot) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot2 <- rbind(output_prot2, result_row)
  }
}

output_prot2$adi_label <- ifelse(output_prot2$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot2$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

plot <- merge(output_prot, output_prot2, by=c("term", "adi_label"))
# plot$sig <- ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Both sig", ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y > 0.05/threshold, "With PCs sig", ifelse(plot$p.value.x > 0.05/threshold & plot$p.value.y <= 0.05/threshold, "No PCs sig", "NS")))

ggplot(plot, aes(x = estimate.x, y = estimate.y, color=adi_label)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black", lty=2) +   # 45 degree 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size=2, alpha=1, show.legend=TRUE) +
  geom_smooth(method = "lm", fullrange = FALSE, alpha=0.2,  color="black", linewidth=0.5) +
  labs(x="Beta with PCs", y="Beta w/o PCs", color="Adiposity Trait") +
  theme_cowplot() + 
  stat_cor(method = "pearson")


#### 4 - enroll_age vs. age_instance2 ####
## Enroll age
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = baseline_prot) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}

output_prot$adi_label <- ifelse(output_prot$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

## MRI age
output_prot2 <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + age_instance2 + mriNum"))
    model <- lm(lm_formula, data = baseline_prot) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot2 <- rbind(output_prot2, result_row)
  }
}

output_prot2$adi_label <- ifelse(output_prot2$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_prot2$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

plot <- merge(output_prot, output_prot2, by=c("term", "adi_label"))
# plot$sig <- ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Both sig", ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y > 0.05/threshold, "With PCs sig", ifelse(plot$p.value.x > 0.05/threshold & plot$p.value.y <= 0.05/threshold, "No PCs sig", "NS")))

ggplot(plot, aes(x = estimate.x, y = estimate.y, color=adi_label)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black", lty=2) +   # 45 degree 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size=2, alpha=1, show.legend=TRUE) +
  geom_smooth(method = "lm", fullrange = FALSE, alpha=0.2,  color="black", linewidth=0.5) +
  labs(x="Beta w/ age_enroll", y="Beta w/ age_mri", color="Adiposity Trait") +
  theme_cowplot() + 
  stat_cor(method = "pearson")


#### 5 - time_between ####
## w/ time_between
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_prot) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}

output_prot$adi_label <- ifelse(output_prot$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

## w/o time_between
output_prot2 <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = baseline_prot) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot2 <- rbind(output_prot2, result_row)
  }
}

output_prot2$adi_label <- ifelse(output_prot2$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_prot2$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

plot <- merge(output_prot, output_prot2, by=c("term", "adi_label"))
# plot$sig <- ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Both sig", ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y > 0.05/threshold, "With PCs sig", ifelse(plot$p.value.x > 0.05/threshold & plot$p.value.y <= 0.05/threshold, "No PCs sig", "NS")))

ggplot(plot, aes(x = estimate.x, y = estimate.y, color=adi_label)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black", lty=2) +   # 45 degree 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size=2, alpha=1, show.legend=TRUE) +
  geom_smooth(method = "lm", fullrange = FALSE, alpha=0.2,  color="black", linewidth=0.5) +
  labs(x="Beta w/ time_between", y="Beta w/o time_between", color="Adiposity Trait", title="Proteomics") +
  theme_cowplot() + facet_wrap(~adi_label) +
  stat_cor(method = "pearson")


baseline_metab <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.n9012.20231201.tsv.gz') # N = 9012
analyte <- colnames(baseline_metab)[which(colnames(baseline_metab) == "Clinical_LDL_C"):which(colnames(baseline_metab) == "Omega_6_pct_PUFA")]

output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_metab) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}

output_prot$adi_label <- ifelse(output_prot$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

## w/o time_between
output_prot2 <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = baseline_metab) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot2 <- rbind(output_prot2, result_row)
  }
}

output_prot2$adi_label <- ifelse(output_prot2$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_prot2$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

plot <- merge(output_prot, output_prot2, by=c("term", "adi_label"))
# plot$sig <- ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Both sig", ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y > 0.05/threshold, "With PCs sig", ifelse(plot$p.value.x > 0.05/threshold & plot$p.value.y <= 0.05/threshold, "No PCs sig", "NS")))

ggplot(plot, aes(x = estimate.x, y = estimate.y, color=adi_label)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black", lty=2) +   # 45 degree 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size=2, alpha=1, show.legend=TRUE) +
  geom_smooth(method = "lm", fullrange = FALSE, alpha=0.2,  color="black", linewidth=0.5) +
  labs(x="Beta w/ time_between", y="Beta w/o time_between", color="Adiposity Trait", title="Metabolomic") +
  theme_cowplot() + facet_wrap(~adi_label) +
  stat_cor(method = "pearson")







#### 6 - time_between stratify  ####
## Enroll age
med_time <- median(baseline_prot$time_between)
a <- subset(baseline_prot, baseline_prot$time_between < med_time)
b <- subset(baseline_prot, baseline_prot$time_between > med_time)

output_prot <- data.frame()
for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = a) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot <- rbind(output_prot, result_row)
  }
}
output_prot$adi_label <- ifelse(output_prot$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_prot$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

output_prot2 <- data.frame()
for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = b) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_prot2 <- rbind(output_prot2, result_row)
  }
}
output_prot2$adi_label <- ifelse(output_prot2$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_prot2$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

plot <- merge(output_prot, output_prot2, by=c("term", "adi_label"))

plot$sig <- ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Both sig", 
                   ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y > 0.05/threshold, "Sig below median", 
                          ifelse(plot$p.value.x > 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Sig above median", "NS")))

ggplot(plot, aes(x = estimate.x, y = estimate.y, color=adi_label)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="#999999", lty=2) +   # 45 degree 
  geom_hline(yintercept = 0, linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = 0, linetype = "dasehd", color = "#999999") +
  geom_point(size=2, alpha=1, show.legend=TRUE) +
  geom_smooth(method = "lm", fullrange = FALSE, alpha=0.2,  color="black", linewidth=0.5) +
  labs(x="Beta below med", y="Beta above med", color="Adiposity Trait", title="Proteomics") +
  theme_cowplot() + 
  stat_cor(method = "pearson")





baseline_metab <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.n9012.20231201.tsv.gz') # N = 9012
analyte <- colnames(baseline_metab)[which(colnames(baseline_metab) == "Clinical_LDL_C"):which(colnames(baseline_metab) == "Omega_6_pct_PUFA")]

med_time <- median(baseline_metab$time_between)
a <- subset(baseline_metab, baseline_metab$time_between < med_time)
b <- subset(baseline_metab, baseline_metab$time_between > med_time)

output_metab <- data.frame()
for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = a) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_metab <- rbind(output_metab, result_row)
  }
}
output_metab$adi_label <- ifelse(output_metab$outcome=="vatadjbmi3", "VATadjBMI", 
                                ifelse(output_metab$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

output_metab2 <- data.frame()
for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + enroll_age + mriNum"))
    model <- lm(lm_formula, data = b) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_metab2 <- rbind(output_metab2, result_row)
  }
}
output_metab2$adi_label <- ifelse(output_metab2$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_metab2$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

plot <- merge(output_metab, output_metab2, by=c("term", "adi_label"))

raw <- fread(file = "/Volumes/medpop_esp2/jdron/projects/cihr_metab/analysis/v1_preNov2023/02_MR/data/raw_jwl_jd.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(raw) = c("id", "outcome", "class", "unit", "lipo_size", "lipo_frac", "lipid_type", "general_type", "lipid_groups", "lipo_groups", "detail_groups","label_name")
raw$outcome = str_replace(raw$outcome, "met-d-IDL_IDL", "met-d-IDL")
raw$outcome = str_replace(raw$outcome, "met-d-", "")
raw <- raw[,-c(1)]

# Exclude ratios and pct
raw <- raw[raw$unit!="ratio",]
met_list = raw$outcome
plot <- merge(plot, raw, by.x="term", by.y="outcome")

rm(raw)

plot$sig <- ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Both sig", 
                   ifelse(plot$p.value.x <= 0.05/threshold & plot$p.value.y > 0.05/threshold, "Sig below median", 
                          ifelse(plot$p.value.x > 0.05/threshold & plot$p.value.y <= 0.05/threshold, "Sig above median", "NS")))

ggplot(plot, aes(x = estimate.x, y = estimate.y, color=sig)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black", lty=2) +   # 45 degree 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size=2, alpha=1, show.legend=TRUE) +
  geom_smooth(method = "lm", fullrange = FALSE, alpha=0.2,  color="black", linewidth=0.5) +
  labs(x="Beta below med", y="Beta above med", color="Adiposity Trait", title="Metabolites") +
  theme_cowplot() + facet_wrap(~adi_label) +
  stat_cor(method = "pearson")
