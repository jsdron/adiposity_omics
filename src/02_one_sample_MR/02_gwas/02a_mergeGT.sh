#!/bin/bash -l
#$ -cwd
#$ -N mergeGT
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 1 -R y -binding linear:1
#$ -l h_vmem=8G
#$ -o log/02a_mergeGT.log
#$ -e log/02a_mergeGT.log
#$ -t 1

source /broad/software/scripts/useuse

### from Buu on Dec. 14, 2022
rm mergelist.txt

for chr in {1..22}; do
  ln -s /broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed ukb_chr${chr}_v2.bed
  ln -s /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim ukb_chr${chr}_v2.bim
  ln -s /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam ukb_chr${chr}_v2.fam
  echo "ukb_chr${chr}_v2" >> mergelist.txt
done

## for chr X
  # ln -s /broad/ukbb/genotype/ukb_cal_chrX_v2.bed ukb_chrX_v2.bed
  # ln -s /broad/ukbb/genotype/ukb_snp_chrX_v2.bim ukb_chrX_v2.bim
  # ln -s /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam ukb_chrX_v2.fam
## I edited the mergelist.txt manually to include chrX
