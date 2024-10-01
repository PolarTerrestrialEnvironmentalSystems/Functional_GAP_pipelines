#!/bin/bash

NR_DB="NR_NCBI"

WORK=${PWD}
OUTDIR="output"
OUT_EGGNOG="out.eggnog"

#INPUT example... See the readme.
INPUT="eggnog_eukaryote_OG_proteins_prokGAP.faa"

#creating your own database using predicted proteins
srun mmseqs createdb ${OUTDIR}/${OUT_EGGNOG}/${INPUT} ${OUTDIR}/${OUT_EGGNOG}/${INPUT}.db

#only eukaryotes searching in the database
srun mmseqs taxonomy ${OUTDIR}/${OUT_EGGNOG}/${INPUT}.db ${INPUT}_taxonomy.result tmp_${INPUT} ${NR_DB} \\
--taxon-list 2759 -e 0.00001 --tax-lineage 1 --lca-ranks species,genus,family,order,class,phylum,kingdom,superkingdom --lca-mode 3 -s 5

#eggnog_eukaryote_OG_proteins_prokGAP
#tabular format output for taxonomy
srun mmseqs createtsv ${OUTDIR}/${OUT_EGGNOG}/${INPUT}.db ${INPUT}_taxonomy.result ${INPUT}_taxonomy.tsv
