#### 0 - LIBRARIES ####

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tableone)
library(stringr)
library(broom)
library(tidyverse)
library(ggtext)
library(ggpubr)
library(gtools)
library(ggrepel)
library(ggbreak)

#### 1 - LOAD DATAFRAMES ####

a <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n4684.20231213.tsv.gz')
a$omic <- "Proteomic"
b <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_metab.n9012.20231213.tsv.gz',
           select=c('term','estimate','std.error','statistic', 'p.value','outcome','adi_label'))
b$omic <- "Metabolomic"

output_all <- rbind(a, b)

threshold = 41+1459 
output_all$sig <- ifelse(output_all$p.value<0.05/threshold, "sig", "ns") 
output_metab$adi_label <- ifelse(output_metab$outcome=="vatadjbmi3", "VATadjBMI", 
                                 ifelse(output_metab$outcome=="asatadjbmi3", "ASATadjBMI", "GFATadjBMI"))

rm(a, b)

#### 2 - SIGNIFICANCE TABLE AND BARPLOT ####

#tatable(output_all$sig)
table(subset(output_all, output_all$omic=="metab")$sig, subset(output_all, output_all$omic=="metab")$outcome)
table(subset(output_all, output_all$omic=="prot")$sig, subset(output_all, output_all$omic=="prot")$outcome)

plot <- output_all %>%
  count(adi_label, omic, sig) %>%
  group_by(adi_label, omic) %>%
  mutate(percentage = n / sum(n)*100)
plot$rounded_per <- round(plot$percentage, 1)

# Counts
ggplot(plot[plot$sig=="sig", ], aes(x = adi_label, y = n, fill = omic)) +
  geom_bar(stat = "identity") +
  labs(x = "Adiposity trait", y = "Significant analytes (n)", fill = "") +
  scale_fill_manual(values=c("#E2A062","#EEA7A8")) +
  theme_cowplot() +
  theme(legend.position = c(0.05, 0.98), legend.justification = c(0, 1)) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) 

# Percentage
ggplot(plot[plot$sig=="sig", ], aes(x = adi_label, y = rounded_per, fill = omic)) +
  geom_bar(stat = "identity", position="dodge") + 
  labs(x = "Adiposity trait", y = "Significant analytes (%)", fill = "") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limit = c(0, 100)) +
  scale_fill_manual(values=c("#E2A062","#EEA7A8")) +
  theme_cowplot() +
  theme(legend.position = c(0.05, 0.98), legend.justification = c(0, 1)) +
  geom_text(aes(label = scales::percent(rounded_per/100)), position = position_dodge(width = 0.9), vjust = -0.5)

#### 3 - MIAMI PLOTS - ALL ####
output_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n4684.20231213.tsv.gz')

link <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/2_olink_protein_map_mr.txt", select=c("Assay", "chr", "gene_start"))
link$Assay <- gsub("-", "_", link$Assay)
link <- distinct(link, Assay, .keep_all = T)
output_prot <- merge(output_prot, link, by.x="term", by.y="Assay", all.x=T)
output_prot$chr <- factor(output_prot$chr, levels=c(1:22, "X"))
output_prot$p_min <- ifelse(output_prot$estimate>0, -log10(output_prot$p.value), log10(output_prot$p.value))

rm(link)

data_cum <- output_prot %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(gene_start)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

output_prot <- output_prot %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = gene_start + bp_add) 

axis_set <- output_prot %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- output_prot %>% 
  filter(p.value == min(p.value)) %>% 
  mutate(ylim = abs(floor(log10(p.value))) + 2) %>% 
  pull(ylim)

sig <- 0.05/threshold

# Loop through this for each adiposity trait

## ASAT
df_sub <- output_prot[output_prot$outcome=="asatadjbmi3",]
df_sub$col <- factor(ifelse(df_sub$p.value<sig & df_sub$term%in%output_prot[output_prot$outcome %in% c("gfatadjbmi3", "vatadjbmi3")&output_prot$p.value<sig,]$term,
                            "A", ifelse(df_sub$p.value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
df_sub$Protein_lab <- ifelse(df_sub$p.value<5e-15, df_sub$term, "")
plot_asat <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                color = col)) +
  geom_hline(yintercept = -log10(sig), color = "#691815", linetype = "solid") + 
  geom_hline(yintercept = log10(sig), color = "#691815", linetype = "solid") + 
  geom_point(alpha = 1, size=1.5) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
  scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#502A6E", "#BB6261")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log<sub>10</sub>(P)") + 
  theme_cowplot() +
  geom_hline(yintercept = 0, color = "black", lty=2) + 
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  ) +
  geom_text_repel(aes(label=Protein_lab), size=2)

## GFAT
df_sub <- output_prot[output_prot$outcome=="gfatadjbmi3",]
df_sub$col <- factor(ifelse(df_sub$p.value<sig & df_sub$term%in%output_prot[output_prot$outcome %in% c("asatadjbmi3", "vatadjbmi3")&output_prot$p.value<sig,]$term,
                            "A", ifelse(df_sub$p.value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
df_sub$Protein_lab <- ifelse(df_sub$p.value<5e-15, df_sub$term, "")
plot_gfat <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                color = col)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "solid") + 
  geom_hline(yintercept = log10(sig), color = "grey40", linetype = "solid") + 
  geom_point(alpha = 1, size=1.5) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
  scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#502A6E", "#3BA4B9")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log<sub>10</sub>(P)") + 
  theme_cowplot() +
  geom_hline(yintercept = 0, color = "black", lty=2) + 
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  ) +
  geom_text_repel(aes(label=Protein_lab), size=2)

## VAT
df_sub <- output_prot[output_prot$outcome=="vatadjbmi3",]
df_sub$col <- factor(ifelse(df_sub$p.value<sig & df_sub$term%in%output_prot[output_prot$outcome %in% c("gfatadjbmi3", "asatadjbmi3")&output_prot$p.value<sig,]$term,
                            "A", ifelse(df_sub$p.value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
df_sub$Protein_lab <- ifelse(df_sub$p.value<5e-15, df_sub$term, "")
plot_vat <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                               color = col)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "solid") + 
  geom_hline(yintercept = log10(sig), color = "grey40", linetype = "solid") + 
  geom_point(alpha = 1, size=1.5) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
  scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#502A6E", "#96C3B1")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log<sub>10</sub>(P)", title="VATadjBMI") + 
  theme_cowplot() +
  geom_hline(yintercept = 0, color = "black", lty=2) + 
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  ) +
  geom_text_repel(aes(label=Protein_lab), size=2)

# Arrange the plots
ggarrange(
  plot_asat + ggtitle("ASATadjBMI"),
  plot_gfat + ggtitle("GFATadjBMI"),
  plot_vat + ggtitle("VATadjBMI"),
  labels=c("A","B","C"),
  ncol = 1, nrow = 3)

#### 4 - MIAMI PLOTS - UNIQUE ####
output_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n4684.20231213.tsv.gz')

link <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/2_olink_protein_map_mr.txt", select=c("Assay", "chr", "gene_start"))
link$Assay <- gsub("-", "_", link$Assay)
link <- distinct(link, Assay, .keep_all = T)
output_prot <- merge(output_prot, link, by.x="term", by.y="Assay", all.x=T)
output_prot$chr <- factor(output_prot$chr, levels=c(1:22, "X"))
output_prot$p_min <- ifelse(output_prot$estimate>0, -log10(output_prot$p.value), log10(output_prot$p.value))

rm(link)

data_cum <- output_prot %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(gene_start)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

output_prot <- output_prot %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = gene_start + bp_add) 

axis_set <- output_prot %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- output_prot %>% 
  filter(p.value == min(p.value)) %>% 
  mutate(ylim = abs(floor(log10(p.value))) + 2) %>% 
  pull(ylim)

threshold = 41+1459 
sig <- 0.05/threshold

df_asat <- subset(output_prot, output_prot$outcome=="asatadjbmi3" & output_prot$p.value<sig)
df_gfat <- subset(output_prot, output_prot$outcome=="gfatadjbmi3" & output_prot$p.value<sig)
df_vat <- subset(output_prot, output_prot$outcome=="vatadjbmi3" & output_prot$p.value<sig)

### ASAT only ###
tmp <- subset(df_asat, !(df_asat$term %in% df_gfat$term) & !(df_asat$term %in% df_vat$term))
a <- ggplot(data=tmp, aes(x=reorder(term, -estimate), y=estimate, color=adi_label)) +
  geom_point(show.legend=FALSE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=FALSE) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261")) 

### GFAT only ###
tmp <- subset(df_gfat, !(df_gfat$term %in% df_asat$term) & !(df_gfat$term %in% df_vat$term))
b <- ggplot(data=tmp, aes(x=reorder(term, -estimate), y=estimate, color=adi_label)) +
  geom_point(show.legend=FALSE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=FALSE) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#3BA4B9")) 

ggarrange(
  a + ggtitle("ASATadjBMI"),
  b + ggtitle("GFATadjBMI"),
  labels=c("A","B"),
  ncol = 2, nrow = 1, widths=c(0.25, 0.75))

### VAT only ###
tmp <- subset(df_vat, !(df_vat$term %in% df_asat$term) & !(df_vat$term %in% df_gfat$term))
tmp <- tmp[order(tmp$estimate, decreasing = TRUE), ]
tmp$split <- rep(seq(1, 8), each = 46)

ggplot(data=tmp, aes(x=reorder(term, -estimate), y=estimate, color=adi_label)) +
  geom_point(show.legend=FALSE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=FALSE) +
  #geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~split, scales="free", ncol=9) +
  coord_flip() +
  scale_color_manual(values=c("#96C3B1")) 

#ggplot(tmp, aes(x = bp_cum, y = estimate, color=adi_label)) +
#  geom_point(alpha = 1, size=1.5) +
#  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
#  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
#                width = 0.2, position = position_dodge(width = 0), show.legend=FALSE) +
#  #scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
#  scale_color_manual(values=c("#96C3B1")) +
#  scale_size_continuous(range = c(0.5,3)) +
#  labs(x = NULL, y = "Effect Estimate (SE)", title="VATadjBMI") + 
#  theme_cowplot() +
#  geom_hline(yintercept = 0, color = "black", lty=2) + 
#  theme( 
#    legend.position = "none",
#    panel.grid.major.x = element_blank(),
#    panel.grid.minor.x = element_blank(),
#    axis.title.y = element_markdown(),
#    axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
#  ) +
#  geom_text_repel(aes(label=term), size=2)

#### 5 - MIAMI PLOTS - SHARED ####
output_prot <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/results/lmResults_adi_prot.n4684.20231213.tsv.gz')

link <- fread("/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/2_olink_protein_map_mr.txt", select=c("Assay", "chr", "gene_start"))
link$Assay <- gsub("-", "_", link$Assay)
link <- distinct(link, Assay, .keep_all = T)
output_prot <- merge(output_prot, link, by.x="term", by.y="Assay", all.x=T)
output_prot$chr <- factor(output_prot$chr, levels=c(1:22, "X"))
output_prot$p_min <- ifelse(output_prot$estimate>0, -log10(output_prot$p.value), log10(output_prot$p.value))

rm(link)

data_cum <- output_prot %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(gene_start)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

output_prot <- output_prot %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = gene_start + bp_add) 

axis_set <- output_prot %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- output_prot %>% 
  filter(p.value == min(p.value)) %>% 
  mutate(ylim = abs(floor(log10(p.value))) + 2) %>% 
  pull(ylim)

threshold = 41+1459 
sig <- 0.05/threshold

df_asat <- subset(output_prot, output_prot$outcome=="asatadjbmi3" & output_prot$p.value<sig)
df_gfat <- subset(output_prot, output_prot$outcome=="gfatadjbmi3" & output_prot$p.value<sig)
df_vat <- subset(output_prot, output_prot$outcome=="vatadjbmi3" & output_prot$p.value<sig)

### SHARED BETWEEN ALL 3 ###
tmp <- subset(output_prot, (output_prot$term %in% df_asat$term) & (output_prot$term %in% df_vat$term) & (output_prot$term %in% df_gfat$term))
ggplot(data=tmp, aes(x=term, y=estimate, color=adi_label)) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=TRUE) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="", title="ASAT & GFAT & VAT") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) 

### ASAT + GFAT ###
tmp <- subset(output_prot, (output_prot$term %in% df_asat$term) & !(output_prot$term %in% df_vat$term) & (output_prot$term %in% df_gfat$term))
ggplot(data=tmp, aes(x=term, y=estimate, color=adi_label, alpha = ifelse(p.value < 0.05/threshold, 1, 0.5))) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=TRUE) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="", title="ASAT & GFAT") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) +
  scale_alpha_continuous(guide = FALSE) # hide alpha legend

### ASAT + VAT ###
tmp <- subset(output_prot, (output_prot$term %in% df_asat$term) & (output_prot$term %in% df_vat$term) & !(output_prot$term %in% df_gfat$term))
ggplot(data=tmp, aes(x=term, y=estimate, color=adi_label, alpha = ifelse(p.value < 0.05/threshold, 1, 0.5))) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=TRUE) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="", title="ASAT & VAT") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) +
  scale_alpha_continuous(guide = FALSE) # hide alpha legend

### GFAT + VAT ###
tmp <- subset(output_prot, !(output_prot$term %in% df_asat$term) & (output_prot$term %in% df_vat$term) & (output_prot$term %in% df_gfat$term))
ggplot(data=tmp, aes(x=term, y=estimate, color=adi_label, alpha = ifelse(p.value < 0.05/threshold, 1, 0.5))) +
  geom_point(show.legend=TRUE) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                width = 0.2, position = position_dodge(width = 0), show.legend=TRUE) +
  geom_hline(yintercept=0, lty=2) + 
  labs(x="", y="Effect Estimate (SE)", color="", title="GFAT & VAT") +
  theme_cowplot() +
  theme(legend.position="top", 
        legend.justification = "center",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#BB6261","#3BA4B9","#96C3B1")) +
  scale_alpha_continuous(guide = FALSE) # hide alpha legend

#### 6 - VOLCANO PLOTS ####

ggplot(output_all[output_all$omic=="Proteomic",], aes(estimate, -log10(p.value), color=sig)) +
  #scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="") +
  labs(y="-log<sub>10</sub>(P)", x="Beta Effect Estimate") +
  geom_point()+
  theme_cowplot() +
  #scale_y_continuous(limits=c(0,10.5)) +
  scale_x_continuous(limits=c(-0.75, 0.75)) +
  guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
  #ggtitle("ASATadjBMI") +
  facet_wrap(~adi_label) +
  theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) #+
  #geom_text_repel(aes(label=ifelse(res[res$outc=="CAD",]$all_sens == "Yes" & (res[res$outc=="CAD",]$pval<0.001|abs(res[res$outc=="CAD",]$b)>0.3), 
                                  # res[res$outc=="CAD",]$exp, "")), size=3, force=10, max.overlaps=100)