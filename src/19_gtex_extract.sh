#!/bin/bash

#########################################
# Script: gtex_extract.sh
# Description: Pipeline to extract eQTLs for the MR candidate proteins.
# Key Outputs:
#   - Summary file of eQTLs for 50 tissues for all MR candidate proteins
#########################################

use Tabix

BED_FILE="/medpop/esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/gene_windows.bed"
GTEX_DIR="/medpop/esp2/projects/GTEx/v10/gtex_v10_eqtl/"
OUTFILE="/medpop/esp2/jdron/projects/adiposity/adiposity_omics/data/gtex/target_protein.all_tissues_eqtls.tsv.gz"

# Write output header
echo -e "chr\tpos\tref\talt\tbuild\tgene_id\tvariant_id\ttss_distance\taf\tma_samples\tma_count\tpval_nominal\tslope\tslope_se\tgene_symbol\ttissue" | gzip -c > "$OUTFILE"

# Loop through GTEx files
for GTEX_FILE in "$GTEX_DIR"/*.tsv.gz; do

  TISSUE=$(basename "$GTEX_FILE" | cut -d'.' -f1)

  while IFS=$'\t' read -r chr start end gene ensembl; do
    tabix "$GTEX_FILE" "$chr:$start-$end" 2>/dev/null | \
    awk -v eid="$ensembl" -v gene="$gene" -v tissue="$TISSUE" '
      BEGIN { FS = OFS = "\t" }
      $1 ~ /^#/ { next }
      {
        split($6, g, ".")
        if (g[1] == eid) {
          print $0, gene, tissue
        }
      }'
  done < "$BED_FILE"
done | gzip -c >> "$OUTFILE"
