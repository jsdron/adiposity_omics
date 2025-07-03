# Libraries


```R
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(cowplot)
library(tidyverse)
# install.packages("forcats")
library(forcats)

library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # this is to prevent .log from being generated after using venn.diagram
```

    
    Attaching package: 'dplyr'
    
    
    The following objects are masked from 'package:data.table':
    
        between, first, last
    
    
    The following objects are masked from 'package:stats':
    
        filter, lag
    
    
    The following objects are masked from 'package:base':
    
        intersect, setdiff, setequal, union
    
    
    -- [1mAttaching packages[22m ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- tidyverse 1.3.2 --
    [32mv[39m [34mtibble [39m 3.2.1     [32mv[39m [34mpurrr  [39m 1.0.1
    [32mv[39m [34mreadr  [39m 2.1.4     [32mv[39m [34mforcats[39m 1.0.0
    -- [1mConflicts[22m -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    [31mx[39m [34mdplyr[39m::[32mbetween()[39m   masks [34mdata.table[39m::between()
    [31mx[39m [34mdplyr[39m::[32mfilter()[39m    masks [34mstats[39m::filter()
    [31mx[39m [34mdplyr[39m::[32mfirst()[39m     masks [34mdata.table[39m::first()
    [31mx[39m [34mdplyr[39m::[32mlag()[39m       masks [34mstats[39m::lag()
    [31mx[39m [34mdplyr[39m::[32mlast()[39m      masks [34mdata.table[39m::last()
    [31mx[39m [34mpurrr[39m::[32mtranspose()[39m masks [34mdata.table[39m::transpose()
    Loading required package: grid
    
    Loading required package: futile.logger
    



    NULL


# Result Files


```R
# Get all the metabolomic variable names in and formatted
raw <- fread(file = "/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/data/raw_jwl_jd.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(raw) = c("id", "outcome", "class", "unit", "lipo_size", "lipo_frac", "lipid_type", "general_type", "lipid_groups", "lipo_groups", "detail_groups","label_name")
raw$outcome = str_replace(raw$outcome, "met-d-IDL_IDL", "met-d-IDL")
raw$outcome = str_replace(raw$outcome, "met-d-", "")
met_list = raw$outcome

# Exclude ratios and pct
raw <- raw[raw$unit!="ratio",]
head(raw,2)
```


<table class="dataframe">
<caption>A data.table: 2 x 12</caption>
<thead>
	<tr><th scope=col>id</th><th scope=col>outcome</th><th scope=col>class</th><th scope=col>unit</th><th scope=col>lipo_size</th><th scope=col>lipo_frac</th><th scope=col>lipid_type</th><th scope=col>general_type</th><th scope=col>lipid_groups</th><th scope=col>lipo_groups</th><th scope=col>detail_groups</th><th scope=col>label_name</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>23460</td><td>Ala</td><td>secondary</td><td>measure</td><td>Ala</td><td>Ala</td><td>Ala</td><td>NL</td><td>Amino_acids</td><td>Metabolome</td><td>Amino acids</td><td>Ala</td></tr>
	<tr><td>23461</td><td>Gln</td><td>secondary</td><td>measure</td><td>Gln</td><td>Gln</td><td>Gln</td><td>NL</td><td>Amino_acids</td><td>Metabolome</td><td>Amino acids</td><td>Gln</td></tr>
</tbody>
</table>



For multiple corrections, typically people correct 0.05 alpha threshold by 41 instead of 249, since a lot of metabolites are highly correlated. This is discussed in this Nightingale tutorial that uses the UKB NMR results ("Statistical significance for multiple testing"): https://nightingalehealth.github.io/ggforestplot/articles/nmr-data-analysis-tutorial.html. 


```R
threshold = 41
```


```R
path = "/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/results/from_jiwoo/reverse/"
```

## Genetic instruments


```R
# IVs <- fread('/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/results/...')
# IVs <- IVs[,-c(1)]
# names(IVs) <- gsub(x = names(IVs), pattern = "met-d-", replacement = "")  
# IVs <- data.frame(IVs)
# IVs <- IVs %>%
#    mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
# head(IVs,2)
```


```R
# # Calculate the number of overlapping elements between each pair of columns

# ### Percentage
# overlap_matrix <- matrix(0, ncol(IVs), ncol(IVs))
# for (i in 1:(ncol(IVs) - 1)) {
#    for (j in (i + 1):ncol(IVs)) {
#        # added the -1 because for some reason it was counting 1 extra? i confirmed in the Excel sheet
#      overlap_matrix[i, j] <- (length(Reduce(intersect, list(IVs[, i], IVs[, j])))-1) / (length(unique(c(IVs[, i], IVs[, j])))-1)  #sum(IVs[, i] %in% IVs[, j])
#      overlap_matrix[j, i] <- overlap_matrix[i, j]
#    }
# }

# for (i in 1:ncol(IVs)) {
#      overlap_matrix[i, i] <- 1
# }

# colnames(overlap_matrix) <- names(IVs)
# rownames(overlap_matrix) <- names(IVs)

# head(overlap_matrix,2)

# options(repr.plot.width = 16, repr.plot.height = 16)
# # Create a heatmap plot
# heatmap(overlap_matrix, col = colorRampPalette(c("white", "blue"))(100),
#         xlab = "", ylab = "", main = "Overlap Between IVs")

# overlap_per <- melt(overlap_matrix)
# colnames(overlap_per)[3] <- c("percent_overlap")

# ### Counts ###
# overlap_matrix_count <- matrix(0, ncol(IVs), ncol(IVs))
# for (i in 1:(ncol(IVs) - 1)) {
#    for (j in (i + 1):ncol(IVs)) {
#      overlap_matrix_count[i, j] <- length(Reduce(intersect, list(IVs[, i], IVs[, j])))-1
#      overlap_matrix_count[j, i] <- overlap_matrix_count[i, j]
#    }
# }

# for (i in 1:ncol(IVs)) {
#      overlap_matrix_count[i, i] <- length(unique(c(IVs[, i])))-1
# }

# colnames(overlap_matrix_count) <- names(IVs)
# rownames(overlap_matrix_count) <- names(IVs)

# counts <- melt(overlap_matrix_count)
# colnames(counts)[3] <- c("unique_overlapping_IVs")

# ### Total ###
# overlap_matrix_count <- matrix(0, ncol(IVs), ncol(IVs))
# for (i in 1:(ncol(IVs) - 1)) {
#    for (j in (i + 1):ncol(IVs)) {
#      overlap_matrix_count[i, j] <- length(unique(c(IVs[, i]))) + length(unique(c(IVs[, j]))) - 2
#      overlap_matrix_count[j, i] <- overlap_matrix_count[i, j]
#    }
# }

# for (i in 1:ncol(IVs)) {
#      overlap_matrix_count[i, i] <- length(unique(c(IVs[, i]))) + length(unique(c(IVs[, i]))) - 2
# }

# colnames(overlap_matrix_count) <- names(IVs)
# rownames(overlap_matrix_count) <- names(IVs)

# total <- melt(overlap_matrix_count)
# colnames(total)[3] <- c("total_combined_IVs")

# ### Total unique ###
# overlap_matrix_count <- matrix(0, ncol(IVs), ncol(IVs))
# for (i in 1:(ncol(IVs) - 1)) {
#    for (j in (i + 1):ncol(IVs)) {
#      overlap_matrix_count[i, j] <- length(union(IVs[, i], IVs[, j]))-1
#      overlap_matrix_count[j, i] <- overlap_matrix_count[i, j]
#    }
# }

# for (i in 1:ncol(IVs)) {
#      overlap_matrix_count[i, i] <- length(unique(c(IVs[, i])))-1
# }

# colnames(overlap_matrix_count) <- names(IVs)
# rownames(overlap_matrix_count) <- names(IVs)

# total_uni <- melt(overlap_matrix_count)
# colnames(total_uni)[3] <- c("total_combined_unique_IVs")


# IV_tot <- merge(counts, total_uni, by=c("Var1", "Var2"))
# IV_tot <- merge(IV_tot, overlap_per, by=c("Var1", "Var2"))
# IV_tot <- merge(IV_tot, total, by=c("Var1", "Var2"))
# head(IV_tot,2)
```


```R
# data_melt <- melt(overlap_matrix)   

# ggp <- ggplot(data_melt, aes(Var1, Var2)) + 
#   geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "black")
# ggp
```

## All MR effects


```R
# Prepare outcomes variable with the names of all the FinnGen outcomes
outcomes = c("CHD","CKD","NAFLD","NASH","T2D","HTN","Obesity")

MR_result_all <- data.frame()

for (i in 1:length(outcomes)){
    outcome <- outcomes[i]
    tmp <- fread(file = paste0(path, outcome, "_mr_reverse_results_20230523.csv"), header = TRUE, stringsAsFactors = FALSE)
    tmp <- tmp[,c(4,2,5:9)]
    colnames(tmp)[2] <- "outcome"
    tmp$outcome = str_replace(tmp$outcome, "met-d-IDL_IDL", "met-d-IDL")
    tmp$outcome = str_replace(tmp$outcome , "met-d-", "")
    tmp$effect <- ifelse(tmp$b>0,"+","-")
    tmp$significant <- tmp$pval<0.05/threshold # bonferroni correction
    tmp <- tmp[tmp$method!="Simple mode",] # Art said that ppl don't usually include all of them
    MR_result_all <- rbind(MR_result_all, tmp)
}

MR_result_all$exposure = str_replace(MR_result_all$exposure, "I9_CHD", "CHD")
MR_result_all$exposure = str_replace(MR_result_all$exposure, "N4_CHRONKIDNEYDIS", "CKD")
MR_result_all$exposure = str_replace(MR_result_all$exposure, "I9_HYPTENS", "HTN")
MR_result_all$exposure = str_replace(MR_result_all$exposure, "E4_OBESITY", "Obesity")

MR_result_all <- subset(MR_result_all, MR_result_all$outcome %in% raw$outcome)
MR_result_all$OR <- exp(MR_result_all$b)
MR_result_all$label <- paste0(MR_result_all$exposure,"_",MR_result_all$outcome)

dim(MR_result_all)
head(MR_result_all)
unique(MR_result_all$outcome)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>3696</li><li>11</li></ol>




<table class="dataframe">
<caption>A data.table: 6 x 11</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>outcome</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>label</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>CHD</td><td>Acetate     </td><td>MR Egger                 </td><td>53</td><td> 0.005881261</td><td>0.03235849</td><td>0.85649686</td><td>+</td><td>FALSE</td><td>1.0058986</td><td>CHD_Acetate     </td></tr>
	<tr><td>CHD</td><td>Acetate     </td><td>Weighted median          </td><td>53</td><td> 0.001749474</td><td>0.01439287</td><td>0.90325426</td><td>+</td><td>FALSE</td><td>1.0017510</td><td>CHD_Acetate     </td></tr>
	<tr><td>CHD</td><td>Acetate     </td><td>Inverse variance weighted</td><td>53</td><td>-0.001033959</td><td>0.01366943</td><td>0.93970528</td><td>-</td><td>FALSE</td><td>0.9989666</td><td>CHD_Acetate     </td></tr>
	<tr><td>CHD</td><td>Acetate     </td><td>Weighted mode            </td><td>53</td><td>-0.006570697</td><td>0.02253242</td><td>0.77174444</td><td>-</td><td>FALSE</td><td>0.9934508</td><td>CHD_Acetate     </td></tr>
	<tr><td>CHD</td><td>Acetoacetate</td><td>MR Egger                 </td><td>53</td><td>-0.081769225</td><td>0.04028699</td><td>0.04762338</td><td>-</td><td>FALSE</td><td>0.9214846</td><td>CHD_Acetoacetate</td></tr>
	<tr><td>CHD</td><td>Acetoacetate</td><td>Weighted median          </td><td>53</td><td>-0.029616008</td><td>0.01403754</td><td>0.03487794</td><td>-</td><td>FALSE</td><td>0.9708182</td><td>CHD_Acetoacetate</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Acetate'</li><li>'Acetoacetate'</li><li>'Acetone'</li><li>'Ala'</li><li>'Albumin'</li><li>'ApoA1'</li><li>'ApoB'</li><li>'bOHbutyrate'</li><li>'Cholines'</li><li>'Citrate'</li><li>'Clinical_LDL_C'</li><li>'Creatinine'</li><li>'DHA'</li><li>'Gln'</li><li>'Glucose'</li><li>'Gly'</li><li>'GlycA'</li><li>'HDL_C'</li><li>'HDL_TG'</li><li>'HDL_PL'</li><li>'HDL_CE'</li><li>'HDL_FC'</li><li>'HDL_L'</li><li>'HDL_P'</li><li>'HDL_size'</li><li>'XL_HDL_P'</li><li>'XL_HDL_L'</li><li>'XL_HDL_PL'</li><li>'XL_HDL_C'</li><li>'XL_HDL_CE'</li><li>'XL_HDL_FC'</li><li>'XL_HDL_TG'</li><li>'L_HDL_P'</li><li>'L_HDL_L'</li><li>'L_HDL_PL'</li><li>'L_HDL_C'</li><li>'L_HDL_CE'</li><li>'L_HDL_FC'</li><li>'L_HDL_TG'</li><li>'M_HDL_P'</li><li>'M_HDL_L'</li><li>'M_HDL_PL'</li><li>'M_HDL_C'</li><li>'M_HDL_CE'</li><li>'M_HDL_FC'</li><li>'M_HDL_TG'</li><li>'S_HDL_P'</li><li>'S_HDL_L'</li><li>'S_HDL_PL'</li><li>'S_HDL_C'</li><li>'S_HDL_CE'</li><li>'S_HDL_FC'</li><li>'S_HDL_TG'</li><li>'His'</li><li>'IDL_P'</li><li>'IDL_L'</li><li>'IDL_PL'</li><li>'IDL_C'</li><li>'IDL_CE'</li><li>'IDL_FC'</li><li>'IDL_TG'</li><li>'Ile'</li><li>'LA'</li><li>'Lactate'</li><li>'LDL_C'</li><li>'LDL_TG'</li><li>'LDL_PL'</li><li>'LDL_CE'</li><li>'LDL_FC'</li><li>'LDL_L'</li><li>'LDL_P'</li><li>'LDL_size'</li><li>'L_LDL_P'</li><li>'L_LDL_L'</li><li>'L_LDL_PL'</li><li>'L_LDL_C'</li><li>'L_LDL_CE'</li><li>'L_LDL_FC'</li><li>'L_LDL_TG'</li><li>'M_LDL_P'</li><li>'M_LDL_L'</li><li>'M_LDL_PL'</li><li>'M_LDL_C'</li><li>'M_LDL_CE'</li><li>'M_LDL_FC'</li><li>'M_LDL_TG'</li><li>'S_LDL_P'</li><li>'S_LDL_L'</li><li>'S_LDL_PL'</li><li>'S_LDL_C'</li><li>'S_LDL_CE'</li><li>'S_LDL_FC'</li><li>'S_LDL_TG'</li><li>'Leu'</li><li>'MUFA'</li><li>'non_HDL_C'</li><li>'Omega_3'</li><li>'Omega_6'</li><li>'Phe'</li><li>'Phosphatidylc'</li><li>'Phosphoglyc'</li><li>'PUFA'</li><li>'Pyruvate'</li><li>'Remnant_C'</li><li>'SFA'</li><li>'Sphingomyelins'</li><li>'Total_C'</li><li>'Total_TG'</li><li>'Total_PL'</li><li>'Total_CE'</li><li>'Total_FC'</li><li>'Total_L'</li><li>'Total_P'</li><li>'Total_FA'</li><li>'Total_BCAA'</li><li>'Tyr'</li><li>'Unsaturation'</li><li>'Val'</li><li>'VLDL_C'</li><li>'VLDL_TG'</li><li>'VLDL_PL'</li><li>'VLDL_CE'</li><li>'VLDL_FC'</li><li>'VLDL_L'</li><li>'VLDL_P'</li><li>'VLDL_size'</li><li>'XXL_VLDL_P'</li><li>'XXL_VLDL_L'</li><li>'XXL_VLDL_PL'</li><li>'XXL_VLDL_C'</li><li>'XXL_VLDL_CE'</li><li>'XXL_VLDL_FC'</li><li>'XXL_VLDL_TG'</li><li>'XL_VLDL_P'</li><li>'XL_VLDL_L'</li><li>'XL_VLDL_PL'</li><li>'XL_VLDL_C'</li><li>'XL_VLDL_CE'</li><li>'XL_VLDL_FC'</li><li>'XL_VLDL_TG'</li><li>'L_VLDL_P'</li><li>'L_VLDL_L'</li><li>'L_VLDL_PL'</li><li>'L_VLDL_C'</li><li>'L_VLDL_CE'</li><li>'L_VLDL_FC'</li><li>'L_VLDL_TG'</li><li>'M_VLDL_P'</li><li>'M_VLDL_L'</li><li>'M_VLDL_PL'</li><li>'M_VLDL_C'</li><li>'M_VLDL_CE'</li><li>'M_VLDL_FC'</li><li>'M_VLDL_TG'</li><li>'S_VLDL_P'</li><li>'S_VLDL_L'</li><li>'S_VLDL_PL'</li><li>'S_VLDL_C'</li><li>'S_VLDL_CE'</li><li>'S_VLDL_FC'</li><li>'S_VLDL_TG'</li><li>'XS_VLDL_P'</li><li>'XS_VLDL_L'</li><li>'XS_VLDL_PL'</li><li>'XS_VLDL_C'</li><li>'XS_VLDL_CE'</li><li>'XS_VLDL_FC'</li><li>'XS_VLDL_TG'</li></ol>



## IVW / Wald results


```R
head(subset(MR_result_all, MR_result_all$exposure=="NASH"))
# NASH only has 1 IV...
```


<table class="dataframe">
<caption>A data.table: 6 x 11</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>outcome</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>label</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NASH</td><td>Acetate     </td><td>Wald ratio</td><td>1</td><td> 0.032746070</td><td>0.008022709</td><td>0.000044700</td><td>+</td><td> TRUE</td><td>1.0332881</td><td>NASH_Acetate     </td></tr>
	<tr><td>NASH</td><td>Acetoacetate</td><td>Wald ratio</td><td>1</td><td>-0.003862483</td><td>0.007961832</td><td>0.627587720</td><td>-</td><td>FALSE</td><td>0.9961450</td><td>NASH_Acetoacetate</td></tr>
	<tr><td>NASH</td><td>Acetone     </td><td>Wald ratio</td><td>1</td><td> 0.000733730</td><td>0.007970090</td><td>0.926650047</td><td>+</td><td>FALSE</td><td>1.0007340</td><td>NASH_Acetone     </td></tr>
	<tr><td>NASH</td><td>Ala         </td><td>Wald ratio</td><td>1</td><td>-0.019976258</td><td>0.007901304</td><td>0.011464156</td><td>-</td><td>FALSE</td><td>0.9802219</td><td>NASH_Ala         </td></tr>
	<tr><td>NASH</td><td>Albumin     </td><td>Wald ratio</td><td>1</td><td> 0.008781185</td><td>0.007972610</td><td>0.270714254</td><td>+</td><td>FALSE</td><td>1.0088199</td><td>NASH_Albumin     </td></tr>
	<tr><td>NASH</td><td>ApoA1       </td><td>Wald ratio</td><td>1</td><td>-0.026896627</td><td>0.007299924</td><td>0.000229145</td><td>-</td><td> TRUE</td><td>0.9734619</td><td>NASH_ApoA1       </td></tr>
</tbody>
</table>




```R
ivw <- subset(MR_result_all, MR_result_all$method=="Inverse variance weighted" | MR_result_all$method=="Wald ratio")
ivw$analysis <- "main"

MR_result <- subset(MR_result_all, (MR_result_all$method=="Inverse variance weighted" & MR_result_all$pval< 0.05/threshold) | (MR_result_all$method=="Wald ratio" & MR_result_all$pval< 0.05/threshold)) # legacy variable name 
MR_result$analysis <- "main"
dim(MR_result)
head(MR_result)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>207</li><li>12</li></ol>




<table class="dataframe">
<caption>A data.table: 6 x 12</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>outcome</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>label</th><th scope=col>analysis</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NAFLD</td><td>Acetate</td><td>Inverse variance weighted</td><td>4</td><td> 0.03427877</td><td>0.010267353</td><td>8.41962e-04</td><td>+</td><td>TRUE</td><td>1.0348731</td><td>NAFLD_Acetate</td><td>main</td></tr>
	<tr><td>NAFLD</td><td>ApoA1  </td><td>Inverse variance weighted</td><td>4</td><td>-0.04384957</td><td>0.008421706</td><td>1.92000e-07</td><td>-</td><td>TRUE</td><td>0.9570979</td><td>NAFLD_ApoA1  </td><td>main</td></tr>
	<tr><td>NAFLD</td><td>HDL_PL </td><td>Inverse variance weighted</td><td>4</td><td>-0.03415058</td><td>0.009623906</td><td>3.87408e-04</td><td>-</td><td>TRUE</td><td>0.9664260</td><td>NAFLD_HDL_PL </td><td>main</td></tr>
	<tr><td>NAFLD</td><td>HDL_L  </td><td>Inverse variance weighted</td><td>4</td><td>-0.03281665</td><td>0.009505818</td><td>5.55890e-04</td><td>-</td><td>TRUE</td><td>0.9677160</td><td>NAFLD_HDL_L  </td><td>main</td></tr>
	<tr><td>NAFLD</td><td>HDL_P  </td><td>Inverse variance weighted</td><td>4</td><td>-0.05730730</td><td>0.014462650</td><td>7.42000e-05</td><td>-</td><td>TRUE</td><td>0.9443038</td><td>NAFLD_HDL_P  </td><td>main</td></tr>
	<tr><td>NAFLD</td><td>M_HDL_P</td><td>Inverse variance weighted</td><td>4</td><td>-0.04367959</td><td>0.007768671</td><td>1.88000e-08</td><td>-</td><td>TRUE</td><td>0.9572606</td><td>NAFLD_M_HDL_P</td><td>main</td></tr>
</tbody>
</table>



## Sensitivity Analysis


```R
sens_res <- subset(MR_result_all, MR_result_all$method!="Inverse variance weighted" & MR_result_all$method!="Wald ratio")
sens_res$analysis <- "sensitivity"
```


```R
# Graph showing all the metabolites, even if they are not significant
options(repr.plot.width=30, repr.plot.height=15)
ggplot(data=MR_result_all, 
       aes(x=outcome, y=b, ymin=(b - 1.96 * se), ymax=(b + 1.96 * se), 
           color=method, alpha = ifelse(pval < 0.05, 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE) + 
        geom_hline(yintercept=0, lty=2) +  
        labs(x="Metabolite Outcome", y="Beta Effect (95% CI)", color="MR Method") + 
        theme_cowplot() +
        theme(legend.position="top", legend.justification = "center") +
        facet_wrap(~exposure, scales = "free_y", drop=TRUE, nrow=length(unique(MR_result_all$exposure))) + # 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_alpha_continuous(guide = FALSE)
```

    Warning message:
    "[1m[22mThe `guide` argument in `scale_*()` cannot be `FALSE`. This was deprecated in ggplot2 3.3.4.
    [36mi[39m Please use "none" instead."



    
![png](output_18_1.png)
    



```R
# Graph showing all the metabolites, even if they are not significant
options(repr.plot.width=30, repr.plot.height=15)
ggplot(data=MR_result_all, 
       aes(x=outcome, y=b, ymin=(b - 1.96 * se), ymax=(b + 1.96 * se), 
           color=method, alpha = ifelse(pval < 0.05/threshold , 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE) + 
        geom_hline(yintercept=0, lty=2) +  
        labs(x="Metabolite Outcome", y="Beta Effect (95% CI)", color="MR Method") + 
        theme_cowplot() +
        theme(legend.position="top", legend.justification = "center") +
        facet_wrap(~exposure, scales = "free_y", drop=TRUE, nrow=length(unique(MR_result_all$exposure))) + # 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_alpha_continuous(guide = FALSE)
```


    
![png](output_19_0.png)
    



```R
# For the primary analysis, we only care about results that are IVW significant and have directionally consistant sens effects
metab_pass <- subset(MR_result_all, MR_result_all$exposure %in% MR_result$exposure)

# Plot showing only metabolites that appear significant in at least one outcome
options(repr.plot.width=25, repr.plot.height=12)
ggplot(data=metab_pass, 
       aes(x=outcome, y=b, ymin=(b - 1.96 * se), ymax=(b + 1.96 * se), 
           color=method, alpha = ifelse(pval < 0.05/threshold, 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE) + 
        geom_hline(yintercept=0, lty=2) +  
        labs(x="Metabolite Outcome", y="Beta Effect (95% CI)", color="MR Method") + 
        theme_cowplot() +
        theme(legend.position="top", legend.justification = "center") +
        facet_wrap(~exposure, scales = "free_y", drop=TRUE, nrow=length(unique(MR_result_all$exposure))) +  
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_alpha_continuous(guide = FALSE)
#         coord_flip() + 
#         theme(legend.position="top", legend.justification = "center") +
#         facet_wrap(~outcome, scales = "free_x", drop=TRUE, ncol=length(unique(plot$outcome))) + # 
#         theme(panel.background = element_rect(fill = "gray99", color = NA)) + 
#         scale_alpha_continuous(guide = FALSE)

```


    
![png](output_20_0.png)
    



```R
# Arrange the data frame by the variable of interest and exposure
merged <- merge(MR_result_all, raw, by="outcome")

plot <- merged %>%
  arrange(detail_groups, desc(outcome)) %>%
  mutate(outcome = factor(outcome, levels = unique(outcome)))

# Then reorder exposure based on its order in the data frame
plot$outcome <- factor(plot$outcome, levels = unique(plot$outcome))

options(repr.plot.width=25, repr.plot.height=12)
ggplot(data=plot[plot$method=="Inverse variance weighted" | plot$method=="Wald ratio",], 
       aes(x=outcome, y=OR, ymin=exp(b - 1.96 * se), ymax=exp(b + 1.96 * se), 
           color=detail_groups, alpha = ifelse(pval < 0.05/threshold, 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE) + 
        geom_hline(yintercept=1, lty=2) +  
        labs(x="Metabolite Outcome", y="Odds Ratio (95% CI)", color="") + 
        theme_cowplot() +
        theme(legend.position="top", legend.justification = "center") +
        facet_wrap(~exposure, scales = "free_y", drop=TRUE, nrow=length(unique(plot$exposure))) +  
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_alpha_continuous(guide = FALSE)
```


    
![png](output_21_0.png)
    



```R
unique(metab_pass$exposure)
# since there isnt anything for CHD, CKD, NASH, trying to plot for them below will fail
outcomes
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'NAFLD'</li><li>'NASH'</li><li>'T2D'</li><li>'HTN'</li><li>'Obesity'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CHD'</li><li>'CKD'</li><li>'NAFLD'</li><li>'NASH'</li><li>'T2D'</li><li>'HTN'</li><li>'Obesity'</li></ol>




```R
# Create a df that only has IVW significant and sensitivity analyses in the right direction
metab_pass <- subset(MR_result_all, MR_result_all$label %in% MR_result$label)
metab_pass <- metab_pass %>%
  group_by(label) %>%
  mutate(directional = ifelse(length(unique(effect)) == 1, 1, 0))
metab_pass <- subset(metab_pass, metab_pass$directional==1)
metab_pass <- subset(metab_pass, metab_pass$label %in% MR_result$label)
metab_pass$analysis <- ifelse(metab_pass$method=="Inverse variance weighted" | metab_pass$method=="Wald ratio", "main", "sensitivity")

dim(metab_pass)
head(metab_pass,2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>645</li><li>13</li></ol>




<table class="dataframe">
<caption>A grouped_df: 2 x 13</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>outcome</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>label</th><th scope=col>directional</th><th scope=col>analysis</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NAFLD</td><td>Acetate</td><td>MR Egger       </td><td>4</td><td>0.10003363</td><td>0.03378517</td><td>0.097645656</td><td>+</td><td>FALSE</td><td>1.105208</td><td>NAFLD_Acetate</td><td>1</td><td>sensitivity</td></tr>
	<tr><td>NAFLD</td><td>Acetate</td><td>Weighted median</td><td>4</td><td>0.03516735</td><td>0.01013919</td><td>0.000523454</td><td>+</td><td> TRUE</td><td>1.035793</td><td>NAFLD_Acetate</td><td>1</td><td>sensitivity</td></tr>
</tbody>
</table>




```R
plot_function <- function(i) {
  plot <- subset(metab_pass, metab_pass$exposure == outcomes[i])
  plot <- merge(plot, raw, by="outcome")   
  gg <- ggplot(data = plot,
         aes(x = outcome, y = b, ymin = (b - 1.96 * se), ymax = (b + 1.96 * se),
             color = method, alpha = ifelse(pval < 0.05 / threshold, 1, 0.3))) +
            # color = method, alpha = ifelse(pval < 0.05, 1, 0.3))) +
            # color = method, alpha = ifelse(analysis == "main", 1, 0.7))) +
            # color = method)) +
    geom_pointrange(position = position_jitter(width = 0.2, height = 0), show.legend = TRUE) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "", y = "Beta effect (95% CI)", title = outcomes[i]) +
    theme_cowplot() +
    theme(legend.position = "top", legend.justification = "center") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_alpha_continuous(guide = FALSE)
    return(gg)
}

options(repr.plot.width = 5, repr.plot.height = 5)
print(plot_function(i = 3))

options(repr.plot.width = 10, repr.plot.height = 5)
print(plot_function(i = 4))

options(repr.plot.width = 15, repr.plot.height = 5)
print(plot_function(i = 5))

options(repr.plot.width = 16, repr.plot.height = 5)
print(plot_function(i = 6))

options(repr.plot.width = 10, repr.plot.height = 5)
print(plot_function(i = 7))
```


    
![png](output_24_0.png)
    



    
![png](output_24_1.png)
    



    
![png](output_24_2.png)
    



    
![png](output_24_3.png)
    



    
![png](output_24_4.png)
    



```R
head(metab_pass)
```


<table class="dataframe">
<caption>A grouped_df: 6 x 13</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>outcome</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>label</th><th scope=col>directional</th><th scope=col>analysis</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NAFLD</td><td>Acetate</td><td>MR Egger                 </td><td>4</td><td> 0.10003363</td><td>0.033785167</td><td>0.097645656</td><td>+</td><td>FALSE</td><td>1.1052081</td><td>NAFLD_Acetate</td><td>1</td><td>sensitivity</td></tr>
	<tr><td>NAFLD</td><td>Acetate</td><td>Weighted median          </td><td>4</td><td> 0.03516735</td><td>0.010139188</td><td>0.000523454</td><td>+</td><td> TRUE</td><td>1.0357930</td><td>NAFLD_Acetate</td><td>1</td><td>sensitivity</td></tr>
	<tr><td>NAFLD</td><td>Acetate</td><td>Inverse variance weighted</td><td>4</td><td> 0.03427877</td><td>0.010267353</td><td>0.000841962</td><td>+</td><td> TRUE</td><td>1.0348731</td><td>NAFLD_Acetate</td><td>1</td><td>main       </td></tr>
	<tr><td>NAFLD</td><td>Acetate</td><td>Weighted mode            </td><td>4</td><td> 0.04320670</td><td>0.010198644</td><td>0.024073345</td><td>+</td><td>FALSE</td><td>1.0441537</td><td>NAFLD_Acetate</td><td>1</td><td>sensitivity</td></tr>
	<tr><td>NAFLD</td><td>ApoA1  </td><td>MR Egger                 </td><td>4</td><td>-0.02911102</td><td>0.039634968</td><td>0.539098198</td><td>-</td><td>FALSE</td><td>0.9713086</td><td>NAFLD_ApoA1  </td><td>1</td><td>sensitivity</td></tr>
	<tr><td>NAFLD</td><td>ApoA1  </td><td>Weighted median          </td><td>4</td><td>-0.04113281</td><td>0.009231469</td><td>0.000008360</td><td>-</td><td> TRUE</td><td>0.9597017</td><td>NAFLD_ApoA1  </td><td>1</td><td>sensitivity</td></tr>
</tbody>
</table>




```R
# Dataframes with nominal significance, robust metabolites for each outcome
# ivw_nom <- subset(MR_result_all, MR_result_all$method=="Inverse variance weighted" & MR_result_all$pval<0.05)

# write.table(MR_result, file="/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/results/MRresults_IVWpass_reverse_2023-06-28.tsv", append = FALSE, sep = "\t", dec = ".", 
            # col.names = TRUE, row.names = FALSE)

```

# Venn Diagrams


```R
# Dataframes with statistically significant, robust metabolites for each outcome
ivw <- subset(metab_pass, metab_pass$method=="Inverse variance weighted")

CHD <- ivw[ivw$exposure=="CHD",]
CKD <- ivw[ivw$exposure=="CKD",]
NAFLD <- ivw[ivw$exposure=="NAFLD",]
NASH <- ivw[ivw$exposure=="NASH",]
T2D <- ivw[ivw$exposure=="T2D",]
HTN <- ivw[ivw$exposure=="HTN",]
Obesity <- ivw[ivw$exposure=="Obesity",]
```


```R
input <- list("Obesity" = Obesity$outcome, 
              "HTN" = HTN$outcome, 
              "T2D" = T2D$outcome, 
              "NAFLD" = NAFLD$outcome) 
venn.plot <- venn.diagram(input, filename = NULL)
options(repr.plot.width=8, repr.plot.height=8)
grid.draw(venn.plot)

overlap <- calculate.overlap(x = input)
overlap
```


<dl>
	<dt>$a6</dt>
		<dd></dd>
	<dt>$a12</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'IDL_L'</li><li>'IDL_PL'</li><li>'M_VLDL_CE'</li><li>'XS_VLDL_C'</li><li>'XS_VLDL_CE'</li></ol>
</dd>
	<dt>$a11</dt>
		<dd></dd>
	<dt>$a5</dt>
		<dd></dd>
	<dt>$a7</dt>
		<dd></dd>
	<dt>$a15</dt>
		<dd></dd>
	<dt>$a4</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Glucose'</li><li>'XL_HDL_P'</li><li>'XL_HDL_L'</li><li>'XL_HDL_PL'</li><li>'XL_HDL_C'</li><li>'XL_HDL_CE'</li><li>'XL_HDL_FC'</li><li>'L_HDL_P'</li><li>'L_HDL_C'</li><li>'L_HDL_CE'</li><li>'L_HDL_FC'</li><li>'Ile'</li><li>'Phe'</li><li>'Total_BCAA'</li><li>'Tyr'</li><li>'Val'</li></ol>
</dd>
	<dt>$a10</dt>
		<dd></dd>
	<dt>$a13</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Clinical_LDL_C'</li><li>'IDL_P'</li><li>'IDL_C'</li><li>'IDL_CE'</li><li>'IDL_FC'</li><li>'LDL_FC'</li><li>'L_LDL_C'</li><li>'L_LDL_FC'</li><li>'Sphingomyelins'</li><li>'Total_C'</li><li>'Total_CE'</li><li>'Total_FC'</li></ol>
</dd>
	<dt>$a8</dt>
		<dd></dd>
	<dt>$a2</dt>
		<dd></dd>
	<dt>$a9</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Albumin'</li><li>'Gly'</li><li>'GlycA'</li><li>'HDL_FC'</li><li>'HDL_size'</li><li>'L_HDL_L'</li><li>'L_HDL_PL'</li></ol>
</dd>
	<dt>$a14</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'ApoB'</li><li>'Cholines'</li><li>'LA'</li><li>'LDL_C'</li><li>'LDL_PL'</li><li>'LDL_CE'</li><li>'LDL_L'</li><li>'LDL_P'</li><li>'L_LDL_P'</li><li>'L_LDL_L'</li><li>'L_LDL_PL'</li><li>'L_LDL_CE'</li><li>'M_LDL_P'</li><li>'M_LDL_L'</li><li>'M_LDL_PL'</li><li>'M_LDL_C'</li><li>'M_LDL_CE'</li><li>'M_LDL_FC'</li><li>'S_LDL_P'</li><li>'S_LDL_L'</li><li>'S_LDL_PL'</li><li>'S_LDL_C'</li><li>'S_LDL_CE'</li><li>'S_LDL_FC'</li><li>'non_HDL_C'</li><li>'Omega_6'</li><li>'Phosphatidylc'</li><li>'Remnant_C'</li><li>'Total_PL'</li><li>'Total_L'</li><li>'VLDL_C'</li><li>'VLDL_CE'</li><li>'M_VLDL_P'</li><li>'M_VLDL_PL'</li><li>'M_VLDL_C'</li><li>'M_VLDL_FC'</li><li>'S_VLDL_PL'</li><li>'S_VLDL_C'</li><li>'S_VLDL_CE'</li><li>'S_VLDL_FC'</li><li>'XS_VLDL_P'</li><li>'XS_VLDL_L'</li><li>'XS_VLDL_PL'</li><li>'XS_VLDL_FC'</li></ol>
</dd>
	<dt>$a1</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Ala'</li><li>'HDL_C'</li><li>'HDL_CE'</li><li>'M_HDL_TG'</li><li>'S_HDL_TG'</li><li>'LDL_size'</li><li>'Leu'</li><li>'Unsaturation'</li><li>'VLDL_TG'</li><li>'VLDL_size'</li><li>'XXL_VLDL_P'</li><li>'XXL_VLDL_L'</li><li>'XXL_VLDL_PL'</li><li>'XXL_VLDL_C'</li><li>'XXL_VLDL_FC'</li><li>'XXL_VLDL_TG'</li><li>'XL_VLDL_P'</li><li>'XL_VLDL_L'</li><li>'XL_VLDL_TG'</li></ol>
</dd>
	<dt>$a3</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Acetate'</li><li>'ApoA1'</li><li>'HDL_P'</li><li>'M_HDL_P'</li><li>'M_HDL_L'</li><li>'M_HDL_PL'</li><li>'M_HDL_C'</li><li>'M_HDL_CE'</li><li>'M_HDL_FC'</li></ol>
</dd>
</dl>




    
![png](output_29_1.png)
    


# Profile Plots


```R
ivw_all <- subset(MR_result_all, MR_result_all$method=="Inverse variance weighted" | MR_result_all$method=="Wald ratio")

MR_result <- merge(ivw_all, raw, by="outcome")
MR_result$pass <- ifelse(MR_result$label %in% metab_pass$label, 1, 0)

forMR <- fread('/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/results/MRresults_IVWpass_foward_2023-06-28.tsv')

head(MR_result,2)
MR_result$label_match <- paste0(MR_result$outcome,"_",MR_result$exposure)
head(MR_result$label)

MR_result$forMR <- ifelse(MR_result$label_match %in% forMR$label, "Bi-directional", "Uni-directional")
table(MR_result$forMR)
```


<table class="dataframe">
<caption>A data.table: 2 x 23</caption>
<thead>
	<tr><th scope=col>outcome</th><th scope=col>exposure</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>...</th><th scope=col>unit</th><th scope=col>lipo_size</th><th scope=col>lipo_frac</th><th scope=col>lipid_type</th><th scope=col>general_type</th><th scope=col>lipid_groups</th><th scope=col>lipo_groups</th><th scope=col>detail_groups</th><th scope=col>label_name</th><th scope=col>pass</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>...</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Acetate</td><td>CHD</td><td>Inverse variance weighted</td><td>53</td><td>-0.001033959</td><td>0.01366943</td><td>0.9397053</td><td>-</td><td>FALSE</td><td>0.9989666</td><td>...</td><td>measure</td><td>Acetate</td><td>Acetate</td><td>Acetate</td><td>NL</td><td>Ketone_bodies</td><td>Metabolome</td><td>Ketone bodies</td><td>Acetate</td><td>0</td></tr>
	<tr><td>Acetate</td><td>CKD</td><td>Inverse variance weighted</td><td> 2</td><td> 0.008960617</td><td>0.01813020</td><td>0.6211387</td><td>+</td><td>FALSE</td><td>1.0090009</td><td>...</td><td>measure</td><td>Acetate</td><td>Acetate</td><td>Acetate</td><td>NL</td><td>Ketone_bodies</td><td>Metabolome</td><td>Ketone bodies</td><td>Acetate</td><td>0</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CHD_Acetate'</li><li>'CKD_Acetate'</li><li>'NAFLD_Acetate'</li><li>'NASH_Acetate'</li><li>'T2D_Acetate'</li><li>'HTN_Acetate'</li></ol>




    
     Bi-directional Uni-directional 
                139            1037 



```R
# Arrange the data frame by the variable of interest and exposure
plot <- MR_result %>%
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

# Flip the order of factors so that they appear on the graph the way I want
plot$label_name <- fct_rev(plot$label_name)

# Function for plotting
plot_function <- function(plot_data) {
  ggplot(data = plot_data[plot_data$method == "Inverse variance weighted", ],
         aes(x = label_name, y = b, ymin = (b - se), ymax = (b + se), color = exposure, alpha = ifelse(pass == 1, 1, 0.5))) +
    geom_pointrange(show.legend = TRUE, position = position_dodge(width = 0.3)) +
    labs(x = "", y = "Beta effect (SE)", color = "Cardiometabolic Exposure") +
    theme_cowplot() +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    facet_wrap(~detail_groups, scales = "free", ncol=length(unique(plot_data$detail_groups))) +
    theme(legend.position = "top", legend.justification = "center") +
    scale_alpha_continuous(guide = FALSE)
}

# options(repr.plot.width=20, repr.plot.height=30)
# plot_function(plot)
options(repr.plot.width=18, repr.plot.height=5)
a <- plot_function(plot[plot$detail_groups=="Cholesterol" | plot$detail_groups=="Cholesteryl esters" 
                        | plot$detail_groups=="Free cholesterol" | plot$detail_groups=="Triglycerides"
                        | plot$detail_groups=="Phospholipids" | plot$detail_groups=="Total lipids",])
a

options(repr.plot.width=3, repr.plot.height=5)
b <- plot_function(plot[plot$detail_groups=="Other lipids",])
b

options(repr.plot.width=9, repr.plot.height=5)
c <- plot_function(plot[plot$detail_groups=="Apolipoproteins" | plot$detail_groups=="Particle concentration" 
                        | plot$detail_groups=="Particle size",])
c

```


    
![png](output_32_0.png)
    



    
![png](output_32_1.png)
    



    
![png](output_32_2.png)
    



```R
options(repr.plot.width=20, repr.plot.height=30)
ggplot(data=plot[plot$method=="Inverse variance weighted",], 
              aes(x=label_name, y=b, ymin=(b - 1.96 * se), ymax=(b + 1.96 * se), color=exposure, alpha = ifelse(pass==1, 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE, position=position_dodge(width=0.3)) +
        labs(x="", y="Beta Effect (95% CI)", color="Cardiometabolic Exposure") + 
        theme_cowplot() + 
        geom_hline(yintercept=0, lty=2) + 
        coord_flip() +  # flip coordinates (puts labels on y axis)
        facet_wrap(~detail_groups, scales = "free") +
        theme(legend.position="top", legend.justification = "center") +
        scale_alpha_continuous(guide = FALSE)
```


    
![png](output_33_0.png)
    



```R
options(repr.plot.width=20, repr.plot.height=30)
ggplot(data=plot[plot$method=="Inverse variance weighted",], 
              aes(x=label_name, y=OR, ymin=exp(b - 1.96 * se), ymax=exp(b + 1.96 * se), color=exposure, alpha = ifelse(pass==1, 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE, position=position_dodge(width=0.3)) +
        labs(x="", y="Odds Ratio (95% CI)", color="Cardiometabolic Exposure") + 
        theme_cowplot() + scale_y_log10() +
        geom_hline(yintercept=1, lty=2) + 
        coord_flip() +  # flip coordinates (puts labels on y axis)
        facet_wrap(~detail_groups, scales = "free") +
        theme(legend.position="top", legend.justification = "center") +
        scale_alpha_continuous(guide = FALSE)
```


    
![png](output_34_0.png)
    



```R
head(plot)
```


<table class="dataframe">
<caption>A data.table: 6 x 25</caption>
<thead>
	<tr><th scope=col>outcome</th><th scope=col>exposure</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th><th scope=col>OR</th><th scope=col>...</th><th scope=col>lipo_frac</th><th scope=col>lipid_type</th><th scope=col>general_type</th><th scope=col>lipid_groups</th><th scope=col>lipo_groups</th><th scope=col>detail_groups</th><th scope=col>label_name</th><th scope=col>pass</th><th scope=col>label_match</th><th scope=col>forMR</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>...</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Val</td><td>CHD  </td><td>Inverse variance weighted</td><td> 53</td><td>-0.006004732</td><td>0.013827374</td><td>6.640967e-01</td><td>-</td><td>FALSE</td><td>0.9940133</td><td>...</td><td>Val</td><td>Val</td><td>NL</td><td>Amino_acids</td><td>Branched-chain_amino_acids</td><td>Amino acids</td><td>Val</td><td>0</td><td>Val_CHD  </td><td>Uni-directional</td></tr>
	<tr><td>Val</td><td>CKD  </td><td>Inverse variance weighted</td><td>  2</td><td>-0.026091351</td><td>0.017513201</td><td>1.362741e-01</td><td>-</td><td>FALSE</td><td>0.9742461</td><td>...</td><td>Val</td><td>Val</td><td>NL</td><td>Amino_acids</td><td>Branched-chain_amino_acids</td><td>Amino acids</td><td>Val</td><td>0</td><td>Val_CKD  </td><td>Uni-directional</td></tr>
	<tr><td>Val</td><td>NAFLD</td><td>Inverse variance weighted</td><td>  4</td><td>-0.005755879</td><td>0.009978653</td><td>5.640616e-01</td><td>-</td><td>FALSE</td><td>0.9942607</td><td>...</td><td>Val</td><td>Val</td><td>NL</td><td>Amino_acids</td><td>Branched-chain_amino_acids</td><td>Amino acids</td><td>Val</td><td>0</td><td>Val_NAFLD</td><td>Bi-directional </td></tr>
	<tr><td>Val</td><td>NASH </td><td>Wald ratio               </td><td>  1</td><td>-0.010717808</td><td>0.007750658</td><td>1.667184e-01</td><td>-</td><td>FALSE</td><td>0.9893394</td><td>...</td><td>Val</td><td>Val</td><td>NL</td><td>Amino_acids</td><td>Branched-chain_amino_acids</td><td>Amino acids</td><td>Val</td><td>0</td><td>Val_NASH </td><td>Uni-directional</td></tr>
	<tr><td>Val</td><td>T2D  </td><td>Inverse variance weighted</td><td>156</td><td> 0.074545815</td><td>0.011758081</td><td>2.300000e-10</td><td>+</td><td> TRUE</td><td>1.0773947</td><td>...</td><td>Val</td><td>Val</td><td>NL</td><td>Amino_acids</td><td>Branched-chain_amino_acids</td><td>Amino acids</td><td>Val</td><td>1</td><td>Val_T2D  </td><td>Uni-directional</td></tr>
	<tr><td>Val</td><td>HTN  </td><td>Inverse variance weighted</td><td>179</td><td> 0.021212875</td><td>0.010737431</td><td>4.820002e-02</td><td>+</td><td>FALSE</td><td>1.0214395</td><td>...</td><td>Val</td><td>Val</td><td>NL</td><td>Amino_acids</td><td>Branched-chain_amino_acids</td><td>Amino acids</td><td>Val</td><td>0</td><td>Val_HTN  </td><td>Bi-directional </td></tr>
</tbody>
</table>




```R
head(tmp)
```


<table class="dataframe">
<caption>A data.table: 6 x 9</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>outcome</th><th scope=col>method</th><th scope=col>nsnp</th><th scope=col>b</th><th scope=col>se</th><th scope=col>pval</th><th scope=col>effect</th><th scope=col>significant</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;lgl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>E4_OBESITY</td><td>Acetate     </td><td>MR Egger                 </td><td>30</td><td>-0.01905085</td><td>0.04014029</td><td>0.638748213</td><td>-</td><td>FALSE</td></tr>
	<tr><td>E4_OBESITY</td><td>Acetate     </td><td>Weighted median          </td><td>30</td><td>-0.04069683</td><td>0.01790421</td><td>0.023024313</td><td>-</td><td>FALSE</td></tr>
	<tr><td>E4_OBESITY</td><td>Acetate     </td><td>Inverse variance weighted</td><td>30</td><td>-0.04642635</td><td>0.01453426</td><td>0.001401849</td><td>-</td><td>FALSE</td></tr>
	<tr><td>E4_OBESITY</td><td>Acetate     </td><td>Weighted mode            </td><td>30</td><td>-0.03809111</td><td>0.01957319</td><td>0.061392350</td><td>-</td><td>FALSE</td></tr>
	<tr><td>E4_OBESITY</td><td>Acetoacetate</td><td>MR Egger                 </td><td>30</td><td> 0.02048907</td><td>0.04074048</td><td>0.618956019</td><td>+</td><td>FALSE</td></tr>
	<tr><td>E4_OBESITY</td><td>Acetoacetate</td><td>Weighted median          </td><td>30</td><td> 0.02675499</td><td>0.01786178</td><td>0.134161798</td><td>+</td><td>FALSE</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=20, repr.plot.height=20)

tmp <- subset(plot, plot$method=="Inverse variance weighted" | plot$method=="Wald ratio" & 
              plot$significant==TRUE & plot$exposure=="CHD")
ggplot(data=tmp, 
              aes(x=label_name, y=OR, ymin=exp(b - 1.96 * se), ymax=exp(b + 1.96 * se), 
                  color=exposure, shape = forMR, alpha = ifelse(pass==1, 1, 0.5))) + # , color=outcome
        geom_pointrange(show.legend=TRUE, position=position_dodge(width=0.3)) +
        labs(x="", y="Odds Ratio (95% CI)", color="Cardiometabolic Exposure", shape = "Causation") + 
        theme_cowplot() + scale_shape_manual(values = c(17, 19)) +
        geom_hline(yintercept=1, lty=2) + 
        #scale_color_manual(values=c("#D12D34", "#FF7A00", "#582998", "#5EACD3", "#65BC77", "#F6D73A", "#396BB0")) +
        # coord_flip() +  # flip coordinates (puts labels on y axis)
        facet_wrap(~detail_groups, scales = "free", drop=TRUE) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(legend.position="top", legend.justification = "center") + guides(colour = guide_legend(nrow = 1)) +
        scale_alpha_continuous(guide = FALSE)
```


    
![png](output_37_0.png)
    



```R
options(repr.plot.width=18, repr.plot.height=11)
ggplot(data = plot[plot$detail_groups!="Amino acids" & plot$detail_groups!="Fatty acids" 
                        & plot$detail_groups!="Fluid balance" & plot$detail_groups!="Glycolysis"
                        & plot$detail_groups!="Inflammation" & plot$detail_groups!="Ketone bodies",],
         aes(x = label_name, y = OR, ymin =exp(b - 1.96 * se), ymax =exp(b + 1.96 * se), 
             color = exposure, shape = forMR, alpha = ifelse(pass == 1, 1, 0.5))) +
    geom_pointrange(show.legend = TRUE, position = position_dodge(width = 0.3)) +
    labs(x = "", y = "", color = "Cardiometabolic Exposure", shape = "Causation") +
    theme_cowplot() + scale_shape_manual(values = c(17, 19)) +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_manual(values=c("#D12D34", "#FF7A00", "#582998", "#5EACD3", "#65BC77", "#F6D73A", "#396BB0")) +
    coord_flip() +
    theme(legend.position = "top", legend.justification = "center") + guides(colour = guide_legend(nrow = 1)) +
    facet_wrap(~detail_groups, scales = "free", ncol=6) +
    #theme(legend.position = "none") +
    scale_alpha_continuous(guide = FALSE)

plot$label_name <- str_replace(plot$label_name , "_", " ")

options(repr.plot.width=18, repr.plot.height=2.75)
ggplot(data = plot[plot$detail_groups=="Amino acids" | plot$detail_groups=="Fatty acids" 
                        | plot$detail_groups=="Fluid balance" | plot$detail_groups=="Glycolysis"
                        | plot$detail_groups=="Inflammation" | plot$detail_groups=="Ketone bodies",],
         aes(x = label_name, y = OR, ymin =exp(b - 1.96 * se), ymax =exp(b + 1.96 * se), 
             color = exposure, shape = forMR, alpha = ifelse(pass == 1, 1, 0.5))) +
    geom_pointrange(show.legend = TRUE, position = position_dodge(width = 0.3)) +
    labs(x = "", y = "Odds Ratio (95% CI)", color = "Cardiometabolic Exposure", shape = "Causation") +
    theme_cowplot() + scale_shape_manual(values = c(17, 19)) +
    geom_hline(yintercept = 1, lty = 2) +
    scale_color_manual(values=c("#D12D34", "#FF7A00", "#582998", "#5EACD3", "#65BC77", "#F6D73A", "#396BB0")) +
    coord_flip() +
    theme(legend.position = "top", legend.justification = "center") + guides(colour = guide_legend(nrow = 1)) +
    facet_wrap(~detail_groups, scales = "free", ncol=6) +
    theme(legend.position = "none") +
    scale_alpha_continuous(guide = FALSE)

```


    
![png](output_38_0.png)
    



    
![png](output_38_1.png)
    



```R
metab_pass$label_match <- paste0(metab_pass$outcome,"_",metab_pass$exposure)

metab_pass$forMR <- ifelse(metab_pass$label_match %in% forMR$label, "Bi-directional", "Uni-directional")

tmp <- metab_pass[metab_pass$method=="Inverse variance weighted" | metab_pass$method=="Wald ratio",]

# Get all the metabolomic variable names in and formatted
chord_anno <- fread(file = "/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/data/chord_jd.txt", header = TRUE, stringsAsFactors = FALSE)

# Exclude ratios and pct
chord_anno <- chord_anno[chord_anno$unit!="ratio",]
chord_anno <- chord_anno[,c("exposure", "general_type")]
colnames(chord_anno)[1] <- "outcome"
head(chord_anno,2)
tmp <- merge(tmp, chord_anno, by="outcome")
tmp <- tmp %>% 
  mutate(outcome = factor(outcome, levels = c('exposure','non_HDL_C','Remnant_C',
                                              'Total_C','Total_CE','Total_FC','Total_P','Total_PL','Total_L','Total_TG','VLDL_C','VLDL_CE','VLDL_FC','VLDL_P','VLDL_size','VLDL_PL','VLDL_L','VLDL_TG','XXL_VLDL_C','XXL_VLDL_CE','XXL_VLDL_FC','XXL_VLDL_L','XXL_VLDL_P','XXL_VLDL_PL','XXL_VLDL_TG','XL_VLDL_C','XL_VLDL_CE','XL_VLDL_FC','XL_VLDL_L','XL_VLDL_P','XL_VLDL_PL','XL_VLDL_TG','L_VLDL_C','L_VLDL_CE','L_VLDL_FC','L_VLDL_L','L_VLDL_P','L_VLDL_PL','L_VLDL_TG','M_VLDL_C','M_VLDL_CE','M_VLDL_FC','M_VLDL_L','M_VLDL_P','M_VLDL_PL','M_VLDL_TG','S_VLDL_C','S_VLDL_CE','S_VLDL_FC','S_VLDL_L','S_VLDL_P','S_VLDL_PL','S_VLDL_TG','XS_VLDL_C','XS_VLDL_CE','XS_VLDL_FC','XS_VLDL_L','XS_VLDL_P','XS_VLDL_PL','XS_VLDL_TG','Clinical_LDL_C','LDL_C','LDL_CE','LDL_FC','LDL_P','LDL_size','LDL_PL','LDL_L','LDL_TG','L_LDL_C','L_LDL_CE','L_LDL_FC','L_LDL_L','L_LDL_P','L_LDL_PL','L_LDL_TG','M_LDL_C','M_LDL_CE','M_LDL_FC','M_LDL_L','M_LDL_P','M_LDL_PL','M_LDL_TG','S_LDL_C','S_LDL_CE','S_LDL_FC','S_LDL_L','S_LDL_P','S_LDL_PL','S_LDL_TG','IDL_C','IDL_CE','IDL_FC','IDL_L','IDL_P','IDL_PL','IDL_TG','HDL_C','HDL_CE','HDL_FC','HDL_P','HDL_size','HDL_PL','HDL_L','HDL_TG','XL_HDL_C','XL_HDL_CE','XL_HDL_FC','XL_HDL_L','XL_HDL_P','XL_HDL_PL','XL_HDL_TG','L_HDL_C','L_HDL_CE','L_HDL_FC','L_HDL_L','L_HDL_P','L_HDL_PL','L_HDL_TG','M_HDL_C','M_HDL_CE','M_HDL_FC','M_HDL_L','M_HDL_P','M_HDL_PL','M_HDL_TG','S_HDL_C','S_HDL_CE','S_HDL_FC','S_HDL_L','S_HDL_P','S_HDL_PL','S_HDL_TG','ApoA1','ApoB','Cholines','Phosphatidylc','Phosphoglyc','Sphingomyelins','Phe','Tyr','Ile','Leu','Total_BCAA','Val','Ala','Gln','Gly','His','Albumin','Creatinine','Citrate','Glucose','Lactate','Pyruvate','Acetate','Acetoacetate','Acetone','bOHbutyrate','GlycA','DHA','LA','MUFA','Omega_3','Omega_6','PUFA','SFA','Total_FA','Unsaturation'))) 
```


<table class="dataframe">
<caption>A data.table: 2 x 2</caption>
<thead>
	<tr><th scope=col>outcome</th><th scope=col>general_type</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>non_HDL_C</td><td>Total lipids</td></tr>
	<tr><td>Remnant_C</td><td>Total lipids</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=22, repr.plot.height=5)
ggplot(data = tmp,
         aes(x = outcome, y = OR, ymin = exp(b - 1.96 * se), ymax = exp(b + 1.96 * se),
             color = exposure, shape = forMR)) + #alpha = ifelse(pval < 0.05 / threshold, 1, 0.5)
    geom_pointrange(position = position_jitter(width = 0.2, height = 0), show.legend = TRUE) +
    geom_hline(yintercept = 1, lty = 2) +
    labs(x = "", y = "Odds Ratio (95% CI)", color = "Cardiometabolic Exposure", shape = "Causation") +
    scale_color_manual(values=c("#582998", "#5EACD3", "#65BC77", "#F6D73A", "#396BB0")) +
    theme_cowplot() + scale_shape_manual(values = c(17, 19)) +
    theme(legend.position = "top", legend.justification = "center",
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA)) +
    #facet_wrap(~general_type, scales = "free", nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_alpha_continuous(guide = FALSE)

ggplot(data = tmp,
         aes(x = outcome, y = b, ymin = (b - se), ymax = (b + se),
             color = exposure, shape = forMR)) + # alpha = ifelse(pval < 0.05 / threshold, 1, 0.5)
    geom_pointrange(position = position_jitter(width = 0.2, height = 0), show.legend = TRUE) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "", y = "Beta (SE)", color = "Cardiometabolic Exposure", shape = "Causation") +
    scale_color_manual(values=c("#582998", "#5EACD3", "#65BC77", "#F6D73A", "#396BB0")) +
    theme_cowplot() + scale_shape_manual(values = c(17, 19)) +
    theme(legend.position = "top", legend.justification = "center",
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA)) +
    #facet_wrap(~general_type, scales = "free", nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_alpha_continuous(guide = FALSE)
```


    
![png](output_40_0.png)
    



    
![png](output_40_1.png)
    



```R
pdf("/medpop/esp2/jdron/projects/cihr_metab/analysis/02_MR/figures/poster_reverseMR_07-07-2023.pdf", 
    width = 22.32, height = 5.07)

ggplot(data = tmp,
         aes(x = outcome, y = b, ymin = (b - se), ymax = (b + se),
             color = exposure, shape = forMR)) + # alpha = ifelse(pval < 0.05 / threshold, 1, 0.5)
    geom_pointrange(position = position_jitter(width = 0.2, height = 0), show.legend = TRUE) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "", y = "Beta (SE)", color = "Cardiometabolic Exposure", shape = "Causation") +
    scale_color_manual(values=c("#582998", "#5EACD3", "#65BC77", "#F6D73A", "#396BB0")) +
    theme_cowplot() + scale_shape_manual(values = c(17, 19)) +
    theme(legend.position = "top", legend.justification = "center",
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA)) +
    #facet_wrap(~general_type, scales = "free", nrow=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_alpha_continuous(guide = FALSE)

dev.off()
```


<strong>cairo_pdf:</strong> 2


# FinnGen vs. CARDIoGRAMC4D
Looking at my FinnGen MR results vs Jiwoo's CARDIoGRAMC4D MR results for CAD


```R
check <- fread('/medpop/esp2/jdron/projects/cihr_metab/from_jiwoo/cad_finngen_cardiogram_comparison.csv')
colnames(check)[1] <- "exposure"
check$exposure = str_replace(check$exposure, "met-d-IDL_IDL", "met-d-IDL")
check$exposure = str_replace(check$exposure, "met-d-", "")
check$pval_consistency <- ifelse(check$pval_consistency==1, "Both significant", "Different significance")
check$pval_consistency <- as.factor(check$pval_consistency)
head(check)
```


<table class="dataframe">
<caption>A data.table: 6 x 10</caption>
<thead>
	<tr><th scope=col>exposure</th><th scope=col>b_fg</th><th scope=col>se_fg</th><th scope=col>pval_fg</th><th scope=col>b_cg</th><th scope=col>se_cg</th><th scope=col>pval_cg</th><th scope=col>pval_consistency</th><th scope=col>b_consistency</th><th scope=col>total_consistency</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>ApoB          </td><td> 0.4946782</td><td>0.06389132</td><td>9.75000e-15</td><td> 0.5166457</td><td>0.06966729</td><td>1.208003e-13</td><td>Both significant      </td><td>1</td><td>1</td></tr>
	<tr><td>ApoB_by_ApoA1 </td><td> 0.5416459</td><td>0.05564653</td><td>2.17000e-22</td><td> 0.5330528</td><td>0.06062903</td><td>1.468686e-18</td><td>Both significant      </td><td>1</td><td>1</td></tr>
	<tr><td>Clinical_LDL_C</td><td> 0.5330459</td><td>0.05355276</td><td>2.43000e-23</td><td> 0.6142225</td><td>0.05697765</td><td>4.276030e-27</td><td>Both significant      </td><td>1</td><td>1</td></tr>
	<tr><td>HDL_C         </td><td>-0.1790881</td><td>0.04699037</td><td>1.38313e-04</td><td>-0.1261611</td><td>0.04816013</td><td>8.802876e-03</td><td>Different significance</td><td>1</td><td>0</td></tr>
	<tr><td>HDL_CE        </td><td>-0.1833821</td><td>0.04901822</td><td>1.83216e-04</td><td>-0.1381596</td><td>0.04889053</td><td>4.714832e-03</td><td>Different significance</td><td>1</td><td>0</td></tr>
	<tr><td>HDL_L         </td><td>-0.1746300</td><td>0.05219519</td><td>8.20718e-04</td><td>-0.1085108</td><td>0.05759519</td><td>5.956153e-02</td><td>Different significance</td><td>1</td><td>0</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=8, repr.plot.height=6)
ggplot(check, aes(x = b_fg, y = b_cg)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = b_cg - se_cg, ymax = b_cg + se_cg), width = 0.01, color="lightgray", alpha=0.7) +
  geom_errorbarh(aes(xmin = b_fg - se_fg, xmax = b_fg + se_fg), height = 0.01, color="lightgray", alpha=0.7) +
  geom_point(aes(color=pval_consistency)) +
  labs(x="FinnGen beta", y="CARDIoGRAM-C4D beta", color="Significance") +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black") + # 45 degree line
  theme_cowplot()
```


    
![png](output_44_0.png)
    



```R
library(cowplot)

options(repr.plot.width=8, repr.plot.height=6)
# assuming x daftaframe is reprieve
ggplot(merged_df, aes(x = beta.x, y = beta.y)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = b_cg - se_cg, ymax = b_cg + se_cg), width = 0.01, color="lightgray", alpha=0.7) +
  geom_errorbarh(aes(xmin = b_fg - se_fg, xmax = b_fg + se_fg), height = 0.01, color="lightgray", alpha=0.7) +
  geom_point(aes(color=pval_consistency)) +
  labs(x="repreieve beta", y="CARDIoGRAM-C4D beta", color="Significance") +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour="black") + # 45 degree line
  theme_cowplot()
```


    Error in ggplot(merged_df, aes(x = beta.x, y = beta.y)): could not find function "ggplot"
    Traceback:




```R
assume
```


```R

```
