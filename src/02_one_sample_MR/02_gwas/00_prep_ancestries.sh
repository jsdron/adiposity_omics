mkdir -p ./tmp/

ancestry=/medpop/esp2/skoyama/passing/ukb_kgpprojection/v02/out/ukb.kgp_projected.tsv.gz

for pop in {EUR,AFR,SAS,EAS,AMR}; do

  zcat $ancestry | grep $pop > ./tmp/$pop.sample

done
