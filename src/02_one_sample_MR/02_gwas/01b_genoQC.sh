#!/bin/bash -l
#$ -cwd
#$ -N clean_QC
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 1 -R y -binding linear:1
#$ -l h_vmem=12G
#$ -o log/01b_clean_QC.log
#$ -e log/01b_clean_QC.log
#$ -t 1

i=$SGE_TASK_ID

source /broad/software/scripts/useuse

# good SNPs
cat ./data/tmp.chr* > ./data/ukb.qc.EUR.snplist

cp /broad/hptmp/jdron/tmp/gwas_qc/ukb.tmp.chr*.mindrem.id ./tmp/

rm ./data/tmp.chr*.all.snplist
