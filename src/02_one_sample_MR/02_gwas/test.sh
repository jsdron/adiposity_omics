#!/bin/bash -l
#$ -cwd
#$ -N filter_rerun
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 1 -R y -binding linear:1
#$ -l h_vmem=12G
#$ -o log/filter_rerun.log
#$ -e log/filter_rerun.log
#$ -t 1-23

i=$SGE_TASK_ID

source /broad/software/scripts/useuse
use PLINK2

chr=$(echo $(echo $(seq 1 22) X) | tr " " "\n" | head -n $i | tail -n 1)

bed=/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed
bim=/broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim
fam=/medpop/esp2/projects/UK_Biobank/linkers/app7089/ukb22418_c1_b0_v2_s488175.fam

tmp=/broad/hptmp/jdron/tmp/gwas_qc/ukb.tmp.NOFILT.chr$chr

### Filter based on missingness
/medpop/esp2/jdron/software/plink2 \
  --bed $bed \
  --bim $bim \
  --fam $fam --no-fid --no-parents --no-sex --no-pheno \
  --write-snplist --write-samples --no-id-header \
  --out $tmp
