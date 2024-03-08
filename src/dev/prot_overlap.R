library(data.table)

adi_pro <- fread('/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/baseline_adi_prot.n4684.20231201.tsv.gz')

all_pro <- fread('/Volumes/medpop_esp2/projects/UK_Biobank/baskets/ukb672544/olink_data_pivot.tsv.gz', select=c("eid"))

pro_only <- subset(all_pro, !(all_pro$eid %in% adi_pro$eid))

write.table(pro_only, "/Volumes/medpop_esp2/jdron/projects/adiposity/adiposity_omics/data/3k_prot_noAdi.tsv",
            quote=FALSE, row.names=FALSE, col.names=TRUE)
