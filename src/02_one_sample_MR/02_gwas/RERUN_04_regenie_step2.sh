#!/bin/bash -l
#$ -cwd
#$ -N reg2_rerun
#$ -l h_rt=144:00:00
#$ -l s_rt=144:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=24G
#$ -o log/reg2_rerun.log
#$ -e log/reg2_rerun.log
#$ -t 1-22

i=$SGE_TASK_ID
chr=$(echo $(echo $(seq 1 22) X) | tr " " "\n" | head -n $i | tail -n 1)

output=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/step2/chr${chr}/

#trait=XXL_VLDL_FC,XL_VLDL_P,XL_VLDL_PL,XL_VLDL_CE,XL_VLDL_FC,XL_VLDL_TG,L_VLDL_P,L_VLDL_PL,L_VLDL_CE,L_VLDL_FC,L_VLDL_TG,M_VLDL_P,M_VLDL_PL,M_VLDL_CE,M_VLDL_FC,M_VLDL_TG,S_VLDL_P,S_VLDL_PL,S_VLDL_CE,S_VLDL_FC,S_VLDL_TG,XS_VLDL_P
trait=Total_FC_by_CE

genData=/broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen # Imputed data
sample=/medpop/esp2/projects/UK_Biobank/linkers/app7089/ukb22828_c1_b0_v3_s487207.sample # Sample file must correspond to input BGEN file

keep_sample=/medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/QCsamplelist_2022-09-08.tsv
covar=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_covar_v1_2022-09-14.tsv
phenofile=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-09-14.tsv #/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-08-29.tsv
## pheno has already undergone rINT

step1=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/step1/combined_pred_added.list

sif=/medpop/esp2/projects/software/singularity/regenie/v3.1.3/regenie.v3.1.3.sif

singularity exec --bind /medpop/:/medpop/,/broad/hptmp/:/broad/hptmp,/broad/ukbb/:/broad/ukbb/ $sif \
  regenie \
    --step 2 \
    --bgen $genData \
    --sample $sample \
    --keep $keep_sample \
    --covarFile $covar \
    --phenoFile $phenofile \
    --bsize 1000 \
    --threads 4 \
    --htp UKB \
    --qt \
    --phenoColList $trait \
    --covarColList enroll_age,fasting,PC{1:20},sex_age,age2,sex_age2 \
    --catCovarList sex,knn,cholmed,fishoil,genotyping_array \
    --minINFO 0.4 \
    --minMAC 10 \
    --pred ${step1} \
    --print-pheno \
    --write-samples \
    --out ${output}ukb.metabs.$chr


#   --phenoColList L_VLDL_P \ ### remove this line for the full run!!!
