#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Salmon Gene Catalog Index
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

# Salmon index
SALMON_INDEX="non_redundant_PROKGAP_protein_to_pCDS.fna.index"

mkdir -p ${OUTDIR}/${OUT_SALMON}

# Here check If there is an index and If not, copy the file to SSD for faster processing.

if [[ -n "$(find ${OUTDIR}/${OUT_SALMON}/ -name '*.index')" ]]
then
        echo "${SALMON_INDEX} found in ${OUTDIR}/${OUT_SALMON}"
else
        rsync -ur ${OUTDIR}/${OUT_PRODIGAL}/${OUT_NUCL} /tmp/.
        srun salmon index -t /tmp/${OUT_NUCL} -i /tmp/${SALMON_INDEX} -p 64
        rsync -ur /tmp/${SALMON_INDEX} ${OUTDIR}/${OUT_SALMON}/.
fi
