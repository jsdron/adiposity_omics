#!/bin/bash -l
#$ -cwd
#$ -N tag
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=64G
#$ -o log/02b_tag.log
#$ -e log/02b_tag.log
#$ -t 1

i=$SGE_TASK_ID

source /broad/software/scripts/useuse

rm /broad/hptmp/jdron/tmp/gwas_qc/ukb.tmp.chr* # generated from 01a

## To merge into a single file, I did:
## OLD USING PREVIOUS FILES FROM MESBAH
# /medpop/esp2/jdron/software/plink2 \
#   --bed /broad/hptmp/mesbah/shared/ukbb/genotyped/ukb_v2_GT_chr1_22.bed \
#   --bim /broad/hptmp/mesbah/shared/ukbb/genotyped/ukb_v2_GT_chr1_22.bim \
#   --fam /broad/hptmp/mesbah/shared/ukbb/genotyped/ukb_v2_GT_chr1_22.fam \
#   --pmerge-list /medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/mergeUKB_bfiles.tsv \
#   --make-bed \
#   --out /broad/hptmp/jdron/tmp/gwas_qc/ukb.merged.genotypes

/medpop/esp2/btruong/Tools/plink --merge-list mergelist.txt --make-bed --memory 64000 --out /broad/hptmp/jdron/tmp/ukb_cal

# i then copied these outputs to my esp2 folder...

# /medpop/esp2/jdron/software/plink2 \
#   --bfile /broad/hptmp/jdron/tmp/ukb_cal \
#   --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
#   --mind 0.1 \
#   --write-snplist --write-samples --no-id-header \
#   --out /broad/hptmp/jdron/tmp/qc_pass

# Copy the output files to a place in /medpop/ so they won't get deleted
# cp /broad/hptmp/jdron/tmp/qc_pass.* /medpop/esp2/jdron/ukbb/genotypes/Dec2022/
# cp /broad/hptmp/jdron/tmp/ukb_cal.log /medpop/esp2/jdron/ukbb/genotypes/Dec2022/



