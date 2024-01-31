#!/bin/bash -l
#$ -cwd
#$ -N reg1
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=6G
#$ -o log/03a_reg1.log
#$ -e log/03a_reg1.log
#$ -t 1

## number above should be the number of analytes being processed #1-168

i=$SGE_TASK_ID

output=/broad/hptmp/jdron/adiposity_omics/metab/gwas/results/step1
mkdir -p ${output}

## INCLUDING X CHROMSOME tag=/broad/hptmp/jdron/tmp/gwas_qc/ukb.merged.genotypes
#Mesbah files tag=/broad/hptmp/mesbah/shared/ukbb/genotyped/ukb_v2_GT_chr1_22
#old tag=/broad/hptmp/jdron/ukb/ukb_v2_GT_chr1_22
tag=/medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/data/ukb_cal

keep_sample=/medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/data/metab_gwas.keep_samples.n82680.20240123.txt
fileredSNPs=/medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/data/ukb.qc.EUR.snplist
covar=/medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/data/metab_gwas.COVAR.n96613.20240123.tsv
## it assumes all variables are quantitative unless otherwise indicated
phenofile=/medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/data/metab_gwas.PHENO.n96613.20240123.tsv

sif=/medpop/esp2/projects/software/singularity/regenie/v3.1.3/regenie.v3.1.3.sif

cat /medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/data/metab_gwas.list.txt | head -n $i | tail -n 1 | while read trait; do

  singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp $sif \
      regenie \
        --step 1 \
        --bed ${tag} \
        --extract $fileredSNPs \
        --covarFile $covar \
        --phenoFile $phenofile \
        --keep $keep_sample \
        --qt \
        --phenoColList $trait \
        --covarColList enroll_age,fasting_time,PC{1:10},enroll_age2 \
        --catCovarList sex,cholmed,fishoil,currentsmoker \
        --bsize 1000 \
        --threads 4 \
        --out ${output}$trait

done
