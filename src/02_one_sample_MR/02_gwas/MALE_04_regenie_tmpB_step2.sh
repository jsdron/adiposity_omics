#!/bin/bash -l
#$ -cwd
#$ -N reg2b_male
#$ -l h_rt=144:00:00
#$ -l s_rt=144:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=24G
#$ -o log/reg2b_male.log
#$ -e log/reg2b_male.log
#$ -t 1-22

i=$SGE_TASK_ID
chr=$(echo $(echo $(seq 1 22) X) | tr " " "\n" | head -n $i | tail -n 1)

output=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step2/chr${chr}/
mkdir -p ${output}

trait=S_VLDL_TG_pct,XS_VLDL_CE_pct,XS_VLDL_FC_pct,XS_VLDL_C_pct,XS_VLDL_PL_pct,XS_VLDL_TG_pct,L_LDL_CE_pct,L_LDL_FC_pct,L_LDL_C_pct,L_LDL_PL_pct,L_LDL_TG_pct,M_LDL_CE_pct,M_LDL_FC_pct,M_LDL_C_pct,M_LDL_PL_pct,M_LDL_TG_pct,S_LDL_CE_pct,S_LDL_FC_pct,S_LDL_C_pct,S_LDL_PL_pct,S_LDL_TG_pct,IDL_CE_pct,IDL_FC_pct,IDL_C_pct,IDL_PL_pct,IDL_TG_pct,XL_HDL_CE_pct,XL_HDL_FC_pct,XL_HDL_C_pct,XL_HDL_PL_pct,XL_HDL_TG_pct,L_HDL_CE_pct,L_HDL_FC_pct,L_HDL_C_pct,L_HDL_PL_pct,L_HDL_TG_pct,M_HDL_CE_pct,M_HDL_FC_pct,M_HDL_C_pct,M_HDL_PL_pct,M_HDL_TG_pct,S_HDL_CE_pct,S_HDL_FC_pct,S_HDL_C_pct,S_HDL_PL_pct,S_HDL_TG_pct,Omega_3_pct,Omega_6_pct,LA_pct,MUFA_pct,PUFA_pct,SFA_pct,DHA_pct,Omega_6_by_Omega_3,ApoB_by_ApoA1,PUFA_by_MUFA,TG_by_PG,XXL_VLDL_FC_pct_C,XL_VLDL_FC_pct_C,L_VLDL_FC_pct_C,M_VLDL_FC_pct_C,S_VLDL_FC_pct_C,XS_VLDL_FC_pct_C,L_LDL_FC_pct_C,M_LDL_FC_pct_C,S_LDL_FC_pct_C,IDL_FC_pct_C,XL_HDL_FC_pct_C,L_HDL_FC_pct_C,M_HDL_FC_pct_C,S_HDL_FC_pct_C,XXL_VLDL_CE_pct_C,XL_VLDL_CE_pct_C,L_VLDL_CE_pct_C,M_VLDL_CE_pct_C,S_VLDL_CE_pct_C,XS_VLDL_CE_pct_C,L_LDL_CE_pct_C,M_LDL_CE_pct_C,S_LDL_CE_pct_C,IDL_CE_pct_C,XL_HDL_CE_pct_C,L_HDL_CE_pct_C,M_HDL_CE_pct_C,S_HDL_CE_pct_C,XXL_VLDL_FC_by_CE,XL_VLDL_FC_by_CE,L_VLDL_FC_by_CE,M_VLDL_FC_by_CE,S_VLDL_FC_by_CE,XS_VLDL_FC_by_CE,L_LDL_FC_by_CE,S_LDL_FC_by_CE,IDL_FC_by_CE,XL_HDL_FC_by_CE,L_HDL_FC_by_CE,M_HDL_FC_by_CE,S_HDL_FC_by_CE,VLDL_CE_pct,VLDL_FC_pct,VLDL_C_pct,VLDL_PL_pct,VLDL_TG_pct,LDL_CE_pct,LDL_FC_pct,LDL_C_pct,LDL_PL_pct,LDL_TG_pct,HDL_CE_pct,HDL_FC_pct,HDL_C_pct,HDL_PL_pct,HDL_TG_pct,VLDL_FC_pct_C,LDL_FC_pct_C,HDL_FC_pct_C,VLDL_CE_pct_C,LDL_CE_pct_C,HDL_CE_pct_C,VLDL_FC_by_CE,LDL_FC_by_CE,HDL_FC_by_CE,Total_CE_pct,Total_FC_pct,Total_C_pct,Total_PL_pct,Total_TG_pct,Total_FC_pct_C,Total_CE_pct_C,Total_FC_by_CE,Omega_3_pct_PUFA,Omega_6_pct_PUFA

genData=/broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen # Imputed data
sample=/medpop/esp2/projects/UK_Biobank/linkers/app7089/ukb22828_c1_b0_v3_s487207.sample # Sample file must correspond to input BGEN file

keep_sample=/medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/QCsamplelist_2022-09-08.tsv
covar=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_covar_v1_2022-09-14.tsv
phenofile=/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-09-14.tsv #/medpop/esp2/jdron/ukbb/pheno_file/gwas/metabolite_pheno_v1_2022-08-29.tsv
## pheno has already undergone rINT

step1=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step1/combined_pred_added.list

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
    --sex-specific 'male' \
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
