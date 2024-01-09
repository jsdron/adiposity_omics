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

baseline_metab <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_metab.n9012.20231201.tsv.gz') # N = 9012

#### 2 - SET VARIABLE LISTS ####

analyte <- colnames(baseline_metab)[which(colnames(baseline_metab) == "Clinical_LDL_C"):which(colnames(baseline_metab) == "Omega_6_pct_PUFA")]

adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

#### 3 - FEATURE ASSOCIATION ####
output_metab <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + age_instance2 + mriNum + time_between"))
    model <- lm(lm_formula, data = baseline_metab) # Change family and response_variable as needed
    summary(model)
    tmp <- tidy(model, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
    result_row <- tmp[which(tmp$term == analyte[j]), ]
    result_row$outcome <- adi_traits[i]
    output_metab <- rbind(output_metab, result_row)
  }
}

output_metab$adi_label <- ifelse(output_metab$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_metab$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

rm(tmp, result_row, i, j)

# For grouping purposes, we need to add the "super groups" for each metabolite. _pct values removed
raw <- fread(file = "/Volumes/medpop_esp2/jdron/projects/cihr_metab/analysis/v1_preNov2023/02_MR/data/raw_jwl_jd.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(raw) = c("id", "outcome", "class", "unit", "lipo_size", "lipo_frac", "lipid_type", "general_type", "lipid_groups", "lipo_groups", "detail_groups","label_name")
raw$outcome = str_replace(raw$outcome, "met-d-IDL_IDL", "met-d-IDL")
raw$outcome = str_replace(raw$outcome, "met-d-", "")
raw <- raw[,-c(1)]

# Exclude ratios and pct
raw <- raw[raw$unit!="ratio",]

met_list = raw$outcome

output_metab <- merge(output_metab, raw, by.x="term", by.y="outcome")

rm(raw)

# Save results
write.table(output_metab, # 9012
            file = gzfile("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n9012.20231213.tsv.gz"), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#### 4 - VISULIZATION ####

# For multiple corrections, typically people correct 0.05 alpha threshold by 41 instead of 249, since a lot of metabolites are highly correlated. This is discussed in this Nightingale tutorial that uses the UKB NMR results ("Statistical significance for multiple testing"): https://nightingalehealth.github.io/ggforestplot/articles/nmr-data-analysis-tutorial.html.
threshold = 41+1459 

output_metab <- output_metab %>%
  arrange(detail_groups, desc(outcome)) %>%
  mutate(label_name = factor(label_name, levels = c('Ala','Gln','Gly','His','Ile','Leu','Phe','Tyr','Val','Total_BCAA',
                                                    'ApoA1','ApoB','HDL','IDL','Clinical LDL','LDL','VLDL','nonHDL',
                                                    'Remnant','Total','CE','FC','C','TG','PL','L','P','DHA','LA',
                                                    'MUFA','Omega_3','Omega_6','PUFA','SFA','Unsaturation','Albumin',
                                                    'Creatinine','Citrate','Glucose','Lactate','Pyruvate','GlycA',
                                                    'Acetate','Acetoacetate','Acetone','bOHbutyrate','Cholines',
                                                    'Phosphatidylc','Phosphoglyc','Sphingomyelins'))) %>%
  mutate(detail_groups = factor(detail_groups, levels = c('Cholesterol','Cholesteryl esters','Free cholesterol',
                                                          'Triglycerides','Phospholipids','Total lipids',
                                                          'Other lipids','Apolipoproteins','Particle concentration',
                                                          'Particle size','HDL - XLarge','HDL - Large','HDL - Medium',
                                                          'HDL - Small','IDL','LDL - Large','LDL - Medium','LDL - Small',
                                                          'Chylomicron','VLDL - XLarge','VLDL - Large','VLDL - Medium',
                                                          'VLDL - Small','VLDL - XSmall','Amino acids','Fatty acids',
                                                          'Fluid balance','Glycolysis','Inflammation','Ketone bodies')))

output_metab <- output_metab[order(output_metab$detail_groups), ]
output_metab$detail_groups <- factor(output_metab$detail_groups)
# Ensure 'detail_groups' is a factor with levels
output_metab$detail_groups <- factor(output_metab$detail_groups)
# Set 'term' as a factor with levels based on 'detail_groups'
output_metab$term <- factor(output_metab$term, levels = unique(plot$term)[order(plot$detail_groups)])

## NO GROUPING ##
plot_all <- ggplot(data=output_metab, aes(x=term, y=estimate, color=adi_label,
                                          alpha = ifelse(p.value < 0.05/threshold, 1, 0.5))) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0)) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="", title="All Metabolites") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) +
  scale_alpha_continuous(guide = FALSE) # hide alpha legend
ggsave("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/metab_plot_all_20231213.pdf", plot_all, width = 25, height = 4)
ggsave("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/metab_plot_all_20231213.png", plot_all, width = 25, height = 4, dpi = 300)

# Significant only
plot <- subset(output_metab, output_metab$p.value<0.05/threshold) 

plot_sig <- ggplot(data=plot, aes(x=term, y=estimate, color=adi_label)) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0)) +
  geom_hline(yintercept=0, lty=2) +  
  labs(x="", y="Effect Estimate (SE)", color="", title="Significant Metabolites") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) +
  scale_alpha_continuous(guide = FALSE) # hide alpha legend
ggsave("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/metab_plot_sig_20231213.pdf", plot_sig, width = 25, height = 4)
ggsave("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/metab_plot_sig_20231213.png", plot_sig, width = 25, height = 4, dpi = 300)

## GROUPING ##
plot_detail <- ggplot(data=output_metab, aes(x=term, y=estimate, color=adi_label,
                                             alpha = ifelse(p.value < 0.05/threshold, 1, 0.5))) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0)) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="", title="All Metabolites") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center") +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) +
  facet_wrap(~detail_groups, scales = "free") + coord_flip() +  # flip coordinates (puts labels on y axis)
  scale_alpha_continuous(guide = FALSE) # hide alpha legend
ggsave("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/metab_plot_all-detail_20231213.pdf", plot_detail, width = 15, height = 10)
ggsave("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/metab_plot_all-detail_20231213.png", plot_detail, width = 15, height = 10, dpi = 300)

#### 5 - TABLE ONE ####
baseline_metab$sex_label <- ifelse(baseline_metab$sex==1, "Male", "Female")

myVars <- c("enroll_age","sex_label",'currentsmoker',"fasting_status","BMI","BMI_class","time_between",
            "asatadjbmi","gfatadjbmi","vatadjbmi")
catVars <- c("sex_label",'currentsmoker',"fasting_status","BMI_class")

tab1 <- CreateTableOne(vars=myVars, data=baseline_metab, factorVars=catVars, addOverall=FALSE)
table1 <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) #nonnormal = nonnorm,
table1

#### 6 - ADIPOSITY HISTOGRAMS ####
a <- ggplot(baseline_metab, aes(x = asatadjbmi3)) +
  geom_density(fill = "#BB6261", color = "black") +
  labs(x = "ASATadjBMI", y = "Density") +
  theme_cowplot()

b <- ggplot(baseline_metab, aes(x = gfatadjbmi3)) +
  geom_density(fill = "#3BA4B9", color = "black") +
  labs(x = "GFATadjBMI", y = "Density") +
  theme_cowplot()

c <- ggplot(baseline_metab, aes(x = vatadjbmi3)) +
  geom_density(fill = "#96C3B1", color = "black") +
  labs(x = "VATadjBMI", y = "Density") +
  theme_cowplot()

ggarrange(a, b, c, labels=c("A", "B", "C"), ncol = 1, nrow = 3)
