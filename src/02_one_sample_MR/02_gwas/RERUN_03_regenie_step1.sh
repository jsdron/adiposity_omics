#!/bin/bash -l
#$ -cwd
#$ -N reg1_rerun
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_vmem=12G
#$ -o log/reg1_rerun.log
#$ -e log/reg1_rerun.log
#$ -t 1

i=$SGE_TASK_ID

output=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/step1/

## INCLUDING X CHROMSOME tag=/broad/hptmp/jdron/tmp/gwas_qc/ukb.merged.genotypes
#Mesbah files tag=/broad/hptmp/mesbah/shared/ukbb/genotyped/ukb_v2_GT_chr1_22
#old tag=/broad/hptmp/jdron/ukb/ukb_v2_GT_chr1_22
tag=/broad/hptmp/jdron/tmp/ukb_cal

keep_sample=/medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/QCsamplelist_2022-09-08.tsv
fileredSNPs=/medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/ukb.qc.all.snplist
covar=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_covar_v1_2022-09-14.tsv
## it assumes all variables are quantitative unless otherwise indicated
phenofile=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-09-14.tsv #/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-08-29.tsv
# pheno has already undergone rINT

sif=/medpop/esp2/projects/software/singularity/regenie/v3.1.3/regenie.v3.1.3.sif

## not including M_LDL_FC_by_CE because all instances were removed from QC

# -t 1-22 # trait=$(echo $(echo XXL_VLDL_FC XL_VLDL_P XL_VLDL_PL XL_VLDL_CE XL_VLDL_FC XL_VLDL_TG L_VLDL_P L_VLDL_PL L_VLDL_CE L_VLDL_FC L_VLDL_TG M_VLDL_P M_VLDL_PL M_VLDL_CE M_VLDL_FC M_VLDL_TG S_VLDL_P S_VLDL_PL S_VLDL_CE S_VLDL_FC S_VLDL_TG XS_VLDL_P) | tr " " "\n" | head -n $i | tail -n 1)
trait=Total_FC_by_CE

  singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp $sif \
      regenie \
        --step 1 \
        --bed ${tag} \
        --extract $fileredSNPs \
        --covarFile $covar \
        --phenoFile $phenofile \
        --phenoColList $trait \
        --keep $keep_sample \
        --qt \
        --covarColList enroll_age,fasting,PC{1:20},sex_age,age2,sex_age2 \
        --catCovarList sex,knn,cholmed,fishoil,genotyping_array \
        --bsize 1000 \
        --threads 4 \
        --out ${output}$trait
