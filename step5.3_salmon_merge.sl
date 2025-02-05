#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Merge Salmon Gene Abundance Results
#===========================================================================

WORK=${PWD}
OUTDIR="output"
OUT_TADPOLE="out.tadpole"
OUT_PRODIGAL="out.prodigal"
OUT_METAEUK="out.metaeuk"
OUT_TIARA="out.tiara"
OUT_SALMON="out.salmon"
OUT_SALMON_MERGED="out.salmon_merged_paired"

if [[-f ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${ID}_conf${CONFIDENCE}_contig.kraken ]]
then
rm -r ${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_all_raw_quant.sf
rm -r ${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_gene_quant.raw.count.len
else
 echo "${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_all_raw_quant.sf and ${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_gene_quant.raw.count.len are not found"

srun python3 ${WORK}/sum_up_qc_merged_paired.py ${OUTDIR}/${OUT_SALMON} ${OUTDIR}/${OUT_SALMON_MERGED}

for file in ${OUTDIR}/${OUT_SALMON_MERGED}/*.sf
do
    filename=$(basename "$file")
    ID="${filename%_merged_paired.quant.sf}"

    # Create a new directory for each ID
    mkdir -p ${OUTDIR}/${OUT_SALMON_MERGED}/${ID}

    # Move the .sf file to the new directory with the new name
    mv "$file" ${OUTDIR}/${OUT_SALMON_MERGED}/${ID}/quant.sf

done

module load salmon/1.10.2

srun salmon quantmerge --quants ${OUTDIR}/${OUT_SALMON_MERGED}/* --column numreads --output ${OUTDIR}/${OUTDIR}_all_raw_quant.sf

srun salmon quantmerge --quants ${OUTDIR}/${OUT_SALMON_MERGED}/* --column len --output ${OUTDIR}/${OUTDIR}_gene_quant.raw.count.len

mv ${OUTDIR}/${OUTDIR}_all_raw_quant.sf ${OUTDIR}/${OUT_SALMON_MERGED}/.
mv ${OUTDIR}/${OUTDIR}_gene_quant.raw.count.len ${OUTDIR}/${OUT_SALMON_MERGED}/.

module unload salmon/1.10.2
