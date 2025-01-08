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


# Salmon index
SALMON_INDEX="non_redundant_PROKGAP_protein_to_pCDS.index"

mkdir -p ${OUTDIR}/${OUT_SALMON}

if [[ -n "$(find ${OUTDIR}/${OUT_PRODIGAL}/ -name '*.index')" ]]
then
        echo "${SALMON_INDEX} found in ${OUTDIR}/${OUT_PRODIGAL}"
else
        srun salmon index -t ${OUTDIR}/${OUT_PRODIGAL}/${OUT_NUCL} -i ${OUTDIR}/${OUT_PRODIGAL}/${SALMON_INDEX}
fi

salmon quant -i ${OUTDIR}/${OUT_PRODIGAL}/${SALMON_INDEX} --meta --libType IU -1 ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_R1.fq.gz -2 ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_R2.fq.gz --validateMappings -o ${OUTDIR}/${OUT_SALMON}/${ID}_paired
salmon quant -i ${OUTDIR}/${OUT_PRODIGAL}/${SALMON_INDEX} --meta --libType U -r ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_merged_R2.fq.gz --validateMappings -o ${OUTDIR}/${OUT_SALMON}/${ID}_merged
