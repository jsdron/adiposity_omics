#!/bin/bash -l
#$ -cwd
#$ -N qc_X
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 1 -R y -binding linear:1
#$ -l h_vmem=12G
#$ -o log/01c_qcX.log
#$ -e log/01c_qcX.log
#$ -t 1

source /broad/software/scripts/useuse
use PLINK2

## These files were generated from 01a_genoQC, so it has already exlcuded any individual and genotype that doesn't match our criteria
bed=/broad/ukbb/genotype/ukb_cal_chrX_v2.bed
bim=/broad/ukbb/genotype/ukb_snp_chrX_v2.bim
fam=/medpop/esp2/projects/UK_Biobank/linkers/app7089/ukb22418_c1_b0_v2_s488175.fam

## this will output a list of all the chromosome X variants, which we can then use in REGENIE to exclude
/medpop/esp2/jdron/software/plink2 \
  --bed $bed \
  --bim $bim \
  --fam $fam --no-fid --no-parents --no-sex --no-pheno \
  --chr X \
  --write-snplist --no-id-header \
  --out ./data/chrX
