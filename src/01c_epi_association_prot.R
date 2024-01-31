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

date <- "20240115"

baseline_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.n4684.20231201.tsv.gz') # N = 4684
colnames(baseline_prot) <- gsub("-", "_", colnames(baseline_prot))

#### 2 - SET VARIABLE LISTS ####

analyte <- colnames(baseline_prot)[which(colnames(baseline_prot) == "CLIP2"):which(colnames(baseline_prot) == "SCARB2")]

adi_traits <- c("vatadjbmi3", "asatadjbmi3", "gfatadjbmi3") # these variables have been adjusted for both BMI and height

#### 3 - FEATURE ASSOCIATION ####
output_prot <- data.frame()

for (i in 1:length(adi_traits)) {
  for (j in 1:length(analyte)) {
    lm_formula <- as.formula(paste0(adi_traits[i], " ~ ", analyte[j], "+ sex + age_instance2 + mriNum + time_between"))
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

rm(tmp, result_row, i, j)

# Save results
write.table(output_prot, 
            file = gzfile(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n4684.",date,".tsv.gz")), 
            append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#### 4 - VISULIZATION ####
threshold = 41+1459 # the number of proteins being considered and 41 from metabs

# Significant only
plot <- subset(output_prot, output_prot$p.value<0.05/threshold) 

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
ggsave(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/prot_plot_sig_",date,".pdf"), plot_sig, width = 25, height = 4)
ggsave(paste0("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/figures/prot_plot_sig_",date,".png"), plot_sig, width = 25, height = 4, dpi = 300)


#### 5 - TABLE ONE ####
baseline_prot$sex_label <- ifelse(baseline_prot$sex==1, "Male", "Female")

myVars <- c("enroll_age","sex_label",'currentsmoker',"fasting_status","BMI","BMI_class","time_between",
            "asatadjbmi","gfatadjbmi","vatadjbmi")
catVars <- c("sex_label",'currentsmoker',"fasting_status","BMI_class")

tab1 <- CreateTableOne(vars=myVars, data=baseline_prot, factorVars=catVars, addOverall=FALSE)
table1 <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) #nonnormal = nonnorm,
table1

#### 6 - ADIPOSITY HISTOGRAMS ####
a <- ggplot(baseline_prot, aes(x = asatadjbmi3)) +
  geom_density(fill = "#BB6261", color = "black") +
  labs(x = "ASATadjBMI", y = "Density") +
  theme_cowplot()

b <- ggplot(baseline_prot, aes(x = gfatadjbmi3)) +
  geom_density(fill = "#3BA4B9", color = "black") +
  labs(x = "GFATadjBMI", y = "Density") +
  theme_cowplot()

c <- ggplot(baseline_prot, aes(x = vatadjbmi3)) +
  geom_density(fill = "#96C3B1", color = "black") +
  labs(x = "VATadjBMI", y = "Density") +
  theme_cowplot()

ggarrange(a, b, c, labels=c("A", "B", "C"), ncol = 1, nrow = 3)