### run this from the same directory that has all of Step1 results

# this will create an empty new file over empty the file if it already exists
echo -n "" > /broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step1/combined_pred_added.list

traitlist=""
lineno=0


### This can be used if you want to do a check and make sure the .loco file exists
### MUST BE RUN FROM THE FOLDER WITH THE LOCO FILES!!!!
while read line; do
    if [ -f ${line}_1.loco ]; then
        echo "${line} /broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step1/${line}_1.loco" >> /broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/male/step1/combined_pred_added.list
        let lineno=lineno+1
        if (( ${lineno} == 1 )); then
            traitlist=${line}
        else
            traitlist=${traitlist}",${line}"
        fi
    fi
done < /medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/allmetab_markers_mod.txt

## save this output file elsewhere so that it doesn't get deleted and I don't have to re-run everything again


### This can be used and it will 100% work
# while read line; do
#         echo "${line} /broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/step1/${line}_1.loco" >> /broad/hptmp/jdron/ukb/metab_gwas/ukb-Metab-Dec2022/result/step1/combined_pred_added.list
# done < /medpop/esp2/jdron/scripts/association_studies/gwas/regenie/runs/ukb-Metab-Aug2022/data/allmetab_markers_mod.txt
