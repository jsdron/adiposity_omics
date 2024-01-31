#!/bin/bash -l
#$ -cwd
#$ -N reg2a_male
#$ -l h_rt=144:00:00
#$ -l s_rt=144:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=24G
#$ -o log/reg2a_male.log
#$ -e log/reg2a_male.log
#$ -t 1-22

i=$SGE_TASK_ID
chr=$(echo $(echo $(seq 1 22) X) | tr " " "\n" | head -n $i | tail -n 1)

output=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step2/chr${chr}/
mkdir -p ${output}

trait=Clinical_LDL_C,VLDL_size,LDL_size,HDL_size,Phosphoglyc,Cholines,Phosphatidylc,Sphingomyelins,ApoB,ApoA1,Unsaturation,Omega_3,Omega_6,MUFA,SFA,LA,DHA,Ala,Gln,Gly,His,Ile,Leu,Val,Phe,Tyr,Glucose,Lactate,Pyruvate,Citrate,bOHbutyrate,Acetate,Acetoacetate,Acetone,Creatinine,Albumin,GlycA,XXL_VLDL_P,XXL_VLDL_PL,XXL_VLDL_CE,XXL_VLDL_FC,XXL_VLDL_TG,XL_VLDL_P,XL_VLDL_PL,XL_VLDL_CE,XL_VLDL_FC,XL_VLDL_TG,L_VLDL_P,L_VLDL_PL,L_VLDL_CE,L_VLDL_FC,L_VLDL_TG,M_VLDL_P,M_VLDL_PL,M_VLDL_CE,M_VLDL_FC,M_VLDL_TG,S_VLDL_P,S_VLDL_PL,S_VLDL_CE,S_VLDL_FC,S_VLDL_TG,XS_VLDL_P,XS_VLDL_PL,XS_VLDL_CE,XS_VLDL_FC,XS_VLDL_TG,IDL_P,IDL_PL,IDL_CE,IDL_FC,IDL_TG,L_LDL_P,L_LDL_PL,L_LDL_CE,L_LDL_FC,L_LDL_TG,M_LDL_P,M_LDL_PL,M_LDL_CE,M_LDL_FC,M_LDL_TG,S_LDL_P,S_LDL_PL,S_LDL_CE,S_LDL_FC,S_LDL_TG,XL_HDL_P,XL_HDL_PL,XL_HDL_CE,XL_HDL_FC,XL_HDL_TG,L_HDL_P,L_HDL_PL,L_HDL_CE,L_HDL_FC,L_HDL_TG,M_HDL_P,M_HDL_PL,M_HDL_CE,M_HDL_FC,M_HDL_TG,S_HDL_P,S_HDL_PL,S_HDL_CE,S_HDL_FC,S_HDL_TG,XXL_VLDL_C,XL_VLDL_C,L_VLDL_C,M_VLDL_C,S_VLDL_C,XS_VLDL_C,IDL_C,L_LDL_C,M_LDL_C,S_LDL_C,XL_HDL_C,L_HDL_C,M_HDL_C,S_HDL_C,XXL_VLDL_L,XL_VLDL_L,L_VLDL_L,M_VLDL_L,S_VLDL_L,XS_VLDL_L,IDL_L,L_LDL_L,M_LDL_L,S_LDL_L,XL_HDL_L,L_HDL_L,M_HDL_L,S_HDL_L,VLDL_CE,VLDL_FC,VLDL_C,VLDL_PL,VLDL_TG,VLDL_L,VLDL_P,LDL_CE,LDL_FC,LDL_C,LDL_PL,LDL_TG,LDL_L,LDL_P,HDL_CE,HDL_FC,HDL_C,HDL_PL,HDL_TG,HDL_L,HDL_P,Total_CE,Total_FC,Total_C,Total_PL,Total_TG,Total_L,Total_P,PUFA,Total_FA,Total_BCAA,non_HDL_C,Remnant_C,XXL_VLDL_CE_pct,XXL_VLDL_FC_pct,XXL_VLDL_C_pct,XXL_VLDL_PL_pct,XXL_VLDL_TG_pct,XL_VLDL_CE_pct,XL_VLDL_FC_pct,XL_VLDL_C_pct,XL_VLDL_PL_pct,XL_VLDL_TG_pct,L_VLDL_CE_pct,L_VLDL_FC_pct,L_VLDL_C_pct,L_VLDL_PL_pct,L_VLDL_TG_pct,M_VLDL_CE_pct,M_VLDL_FC_pct,M_VLDL_C_pct,M_VLDL_PL_pct,M_VLDL_TG_pct,S_VLDL_CE_pct,S_VLDL_FC_pct,S_VLDL_C_pct,S_VLDL_PL_pct

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
