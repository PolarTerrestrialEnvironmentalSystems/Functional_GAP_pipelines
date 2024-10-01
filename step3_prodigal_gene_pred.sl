#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Prodigal Gene Prediction
#===========================================================================

WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_TADPOLE="out.tadpole"
OUT_MEGAHIT="out.megahit"
OUT_PRODIGAL="out.prodigal"
TMP="tmp"

CPU=${SLURM_CPUS_PER_TASK}

# Define extension
END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"

cd ${OUTDIR}/${OUT_TADPOLE}
FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
ID=${FILE_MERGED%${END_MERGED}}
cd ${WORK}

mkdir -p ${OUTDIR}/${OUT_PRODIGAL}

prodigal -g 1 -p meta -i ${OUTDIR}/${OUT_MEGAHIT}/${ID}/final.contigs.fa -a ${OUTDIR}/${OUT_PRODIGAL}/${ID}/${ID}.faa -d  ${OUTDIR}/${OUT_PRODIGAL}/${ID}/${ID}.fna
