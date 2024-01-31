#!/bin/bash -l
#$ -cwd
#$ -N merge_rerun
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -pe smp 1 -R y -binding linear:1
#$ -l h_vmem=12G
#$ -o log/merge_rerun.log
#$ -e log/merge_rerun.log
#$ -t 1-22

i=$SGE_TASK_ID

# there are 22 metabolites listed here
trait=$(echo $(echo XXL_VLDL_FC XL_VLDL_P XL_VLDL_PL XL_VLDL_CE XL_VLDL_FC XL_VLDL_TG L_VLDL_P L_VLDL_PL L_VLDL_CE L_VLDL_FC L_VLDL_TG M_VLDL_P M_VLDL_PL M_VLDL_CE M_VLDL_FC M_VLDL_TG S_VLDL_P S_VLDL_PL S_VLDL_CE S_VLDL_FC S_VLDL_TG XS_VLDL_P) | tr " " "\n" | head -n $i | tail -n 1)

mkdir -p ../result/all/${trait}/

# Where are the results

output=../result/all/${trait}

# Merge chromosome outputs together
for chr in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}; do

  results=/broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/step2/chr${chr}/ukb.metabs.${chr}_${trait}.regenie
  cat $results | sed '1d' >> $output/ukb.${trait}.merged.chr1-chrX.regenie
  cat $results.ids | wc -l >> $output/ukb.${trait}.merged.chr1-chrX.regenie.ids

done

cat $output/ukb.${trait}.merged.chr1-chrX.regenie | sed '1s/^/Name\tChr\tPos\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tInfo\n/' | gzip > $output/ukb.${trait}.merged.chr1-chrX.regenie.gz

rm $output/ukb.${trait}.merged.chr1-chrX.regenie
