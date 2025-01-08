#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Salmon Gene Abundance Estimation
#===========================================================================

# this is a prokGAP input example. This should be changed based on the pipeline.
OUT_NUCL="out.prodigal/non_redundant_PROKGAP_protein_to_pCDS.fna"

WORK=${PWD}
OUTDIR="output"
OUT_TADPOLE="out.tadpole"
OUT_PRODIGAL="out.prodigal"
OUT_METAEUK="out.metaeuk"
OUT_TIARA="out.tiara"
OUT_SALMON="out.salmon"

END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"
CPU=${SLURM_CPUS_PER_TASK}

#===================================================================
# environment
#===================================================================
cd ${WORK}/${OUTDIR}/${OUT_TADPOLE}
FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
ID=${FILE_MERGED%${END_MERGED}}
cd ${WORK}

salmon quant -i ${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX} --meta --libType IU -1 ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_R1.fq.gz -2 ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_R2.fq.gz --validateMappings -o ${OUTDIR}/${OUT_SALMON}/${ID}_paired
salmon quant -i ${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX} --meta --libType U -r ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_merged_R2.fq.gz --validateMappings -o ${OUTDIR}/${OUT_SALMON}/${ID}_merged
