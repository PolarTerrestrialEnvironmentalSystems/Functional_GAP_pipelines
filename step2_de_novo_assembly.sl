#!/bin/bash
#===========================================================================
# Author: Lars Harms
# Modified : Ugur Cabuk
# contact: lars.harms@awi.de, ugur.cabuk@awi.de
# slurm options and variables under >set required variables<
# have to be modified by the user
#=============================================================================


# given variables
#===========================================================================
WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_TADPOLE="out.tadpole"
OUT_MEGAHIT="out.megahit"
TMP="tmp"


END_R1="_fastp_R1.fq.gz"
END_R2="_fastp_R2.fq.gz"
END_MERGED="_fastp_merged_R2.fq.gz"

CPU=${SLURM_CPUS_PER_TASK}

# prepare environment
#===========================================================================
mkdir -p ${OUTDIR}/${OUT_TADPOLE}
mkdir -p ${OUTDIR}/${OUT_MEGAHIT}

mkdir -p ${TMP}

cd ${OUTDIR}/${OUT_FASTP}
FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
ID=${FILE_MERGED%${END_MERGED}}
cd ${WORK}

# tasks to be performed
#===========================================================================

# tadpole
# ------------
tadpole.sh in=${OUTDIR}/${OUT_FASTP}/${FILE_MERGED} out=${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_MERGED} mode=correct k=50
tadpole.sh in=${OUTDIR}/${OUT_FASTP}/${ID}${END_R1} in2=${OUTDIR}/${OUT_FASTP}/${ID}${END_R2} out=${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R1} out2=${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R2} mode=correct k=50

# megahit
#------------
echo "Assembly is performed with megahit."
megahit --presets meta-large -m 1 -t ${CPU} --tmp-dir ${TMP} --continue --min-contig-len 300 -1 ${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R1} -2 ${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R2} -r ${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_MERGED} -o ${OUTDIR}/${OUT_MEGAHIT}/${ID}
