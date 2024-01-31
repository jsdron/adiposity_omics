#!/bin/bash -l
#$ -cwd
#$ -N genoQC
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_vmem=12G
#$ -o log/01a_genoQC.log
#$ -e log/01a_genoQC.log
#$ -t 1-23

i=$SGE_TASK_ID

source /broad/software/scripts/useuse
use PLINK2

chr=$(echo $(echo $(seq 1 22) X) | tr " " "\n" | head -n $i | tail -n 1)

bed=/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed
bim=/broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim
fam=/medpop/esp2/projects/UK_Biobank/linkers/app7089/ukb22418_c1_b0_v2_s488175.fam

## Make sure this tmp folder is present
tmp=/broad/hptmp/jdron/tmp/gwas_qc/ukb.tmp.chr$chr

### Filter based on missingness
/medpop/esp2/jdron/software/plink2 \
  --bed $bed \
  --bim $bim \
  --fam $fam --no-fid --no-parents --no-sex --no-pheno \
  --mind 0.2 \
  --geno 0.2 \
  --keep /medpop/esp2/jdron/projects/adiposity/adiposity_omics/src/02_one_sample_MR/02_gwas/tmp/samplelist.EUR.txt \
  --write-snplist --write-samples --no-id-header \
  --out $tmp.delete

### Filter based on missingness and MAF
/medpop/esp2/jdron/software/plink2 \
  --bed $bed \
  --bim $bim \
  --fam $fam --no-fid --no-parents --no-sex --no-pheno \
  --mind 0.02 \
  --geno 0.01 \
  --maf 0.001 \
  --hwe 1e-6 \
  --keep $tmp.delete.id \
  --extract $tmp.delete.snplist \
  --write-snplist --write-samples --no-id-header \
  --out $tmp

rm $tmp.delete.*

cat /broad/hptmp/jdron/tmp/gwas_qc/ukb.tmp.chr$chr.snplist > ./data/tmp.chr$chr.all.snplist
