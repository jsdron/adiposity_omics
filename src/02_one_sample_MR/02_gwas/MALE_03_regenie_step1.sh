#!/bin/bash -l
#$ -cwd
#$ -N reg1_male
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=6G
#$ -o log/reg1_male.log
#$ -e log/reg1_male.log
#$ -t 1-324

i=$SGE_TASK_ID #1-324

output=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step1/
mkdir -p ${output}

tag=/broad/hptmp/jdron/tmp/ukb_cal

keep_sample=/medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/QCsamplelist_2022-09-08.tsv
fileredSNPs=/medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/ukb.qc.all.snplist
covar=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_covar_v1_2022-09-14.tsv
## it assumes all variables are quantitative unless otherwise indicated
phenofile=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-09-14.tsv #/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-08-29.tsv
# pheno has already undergone rINT

sif=/medpop/esp2/projects/software/singularity/regenie/v3.1.3/regenie.v3.1.3.sif

## not including M_LDL_FC_by_CE because all instances were removed from QC
cat /medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/allmetab_markers_mod.txt | head -n $i | tail -n 1 | while read trait; do

# trait=Total_FC_by_CE # used this temporarily with -t 1 in order to run the trait that got missed

  singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp $sif \
      regenie \
        --step 1 \
        --bed ${tag} \
        --extract $fileredSNPs \
        --exclude /medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/chrX.qc.snplist \
        --covarFile $covar \
        --phenoFile $phenofile \
        --keep $keep_sample \
        --qt \
        --sex-specific 'male' \
        --phenoColList $trait \
        --covarColList enroll_age,fasting,PC{1:20},sex_age,age2,sex_age2 \
        --catCovarList sex,knn,cholmed,fishoil,genotyping_array \
        --bsize 1000 \
        --threads 4 \
        --out ${output}$trait

done

### for sex-specific analysis, include --sex-specific ['female'/'male']
# regenie assumes: 0/1/2 [0=unknown, 1=male, 2=female]
