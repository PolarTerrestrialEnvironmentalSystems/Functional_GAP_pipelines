#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Tiara Contig Classification
#===========================================================================

WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_TADPOLE="out.tadpole"
OUT_MEGAHIT="out.megahit"
OUT_TIARA="out.tiara"
TMP="tmp"

CPU=${SLURM_CPUS_PER_TASK}

# Define extension
END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"

cd ${OUTDIR}/${OUT_TADPOLE}
FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
ID=${FILE_MERGED%${END_MERGED}}
cd ${WORK}

mkdir -p ${OUTDIR}/${OUT_TIARA}

tiara -i ${OUTDIR}/${OUT_MEGAHIT}/${ID}/final.contigs.fa -o ${OUTDIR}/${OUT_TIARA/${ID}_300.fna --min_len 300 --to_fasta arc euk bac pro

cat {OUTDIR}/${OUT_TIARA/archaea_${ID}_300.fna {OUTDIR}/${OUT_TIARA/bacteria_${ID}_300.fna {OUTDIR}/${OUT_TIARA/prokarya_${ID}_300.fna > {OUTDIR}/${OUT_TIARA/prokaryotic_{ID}_300.bin.fna
