#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Metaeuk Gene Prediction
#===========================================================================

WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_TADPOLE="out.tadpole"
OUT_MEGAHIT="out.megahit"
OUT_METAEUK="out.metaeuk"
TMP="tmp"

CPU=${SLURM_CPUS_PER_TASK}

# Database
MMSEQS2_DB="uniref50_database"

# Define extension
END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"

cd ${OUTDIR}/${OUT_TADPOLE}
FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
ID=${FILE_MERGED%${END_MERGED}}
cd ${WORK}

mkdir -p ${OUTDIR}/${OUT_METAEUK}

metaeuk easy-predict ${OUTDIR}/${OUT_MEGAHIT}/${ID}/final.contigs.fa ${MMSEQS2_DB2} ${OUTDIR}/${OUT_METAEUK}/${ID} ${OUTDIR}/${TMP}_${ID} --min-ungapped-score 35 --min-exon-aa 20 --min-length 40 --remove-tmp-files 1




