#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Pydamage Workflow.v1
#===========================================================================

# given variables
#===========================================================================
WORK=${PWD}
OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_TADPOLE="out.tadpole"
OUT_MEGAHIT="out.megahit"
OUT_BWA="out.bwa"
OUT_PYDAMAGE="out.pydamage"
TMP="tmp"

# Set "YES" for what modules you want to run. E.g, If you get error in pydamage, do not need run bwa again, so leave it as a blank.
RUN_BWA="YES"
RUN_PYDAMAGE="YES"
RUN_KRAKEN="YES"

mkdir -p ${WORK}/${OUTDIR}/${OUT_BWA}
mkdir -p ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}

END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"
REF="final.contigs.fa"

cd ${OUTDIR}/${OUT_TADPOLE}
FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
SAMPLE_ID=${FILE_MERGED%${END_MERGED}}
cd ${WORK}

# This parameter for HPC, otherwise remove it here and also from kraken and pydamage parameter.
CPU=${SLURM_CPUS_PER_TASK}

if [ "${RUN_BWA}" = "YES" ]; then
# Here If statement has to be written.
#MODULES
module load bwa/0.7.18
module load samtools/1.20
module load bamtools/2.5.2
module load pydamage/0.80

srun bwa index ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/final.contigs.fa

# BWA-mapping to get coverage
srun bwa mem ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/${REF} ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${SAMPLE_ID}_tadpole_ecc_fastp_merged_R2.fq.gz > ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_merged.sam
srun bwa mem ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/${REF} ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${SAMPLE_ID}_tadpole_ecc_fastp_R1.fq.gz ${OUT_TADPOLE}/${SAMPLE_ID}_tadpole_ecc_fastp_R2.fq.gz > ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_paired.sam

srun samtools sort ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_merged.sam -o ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_merged.sorted.bam
srun samtools sort ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_paired.sam -o ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_paired.sorted.bam

# Remove temporary files
rm -rf ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_merged.sam
rm -rf ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_paired.sam

# Combine merge and paired bam files.

srun bamtools merge -in ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_merged.sorted.bam -in ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_paired.sorted.bam -out ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}.merge_paired.bam

# Remove temporary files
rm -rf ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_merged.sorted.bam
rm -rf ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}_out_paired.sorted.bam

# Sort it again merge bam files
srun samtools sort ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}.merge_paired.bam > ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}.merge_paired.sorted.bam

# Remove temporary file
rm -rf ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}.merge_paired.bam
else
    echo "Skipping BWA ANALYSIS."
fi

# PYDAMAGE
if [ "${RUN_PYDAMAGE}" = "YES" ]; then

srun pydamage --outdir ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID} analyze ${WORK}/${OUTDIR}/${OUT_BWA}/${SAMPLE_ID}.merge_paired.sorted.bam -p ${CPU}

awk -F, -v OFS=, -v prefix="$SAMPLE_ID" 'NR==1{print; next} {print prefix "," $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36}' "${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}/pydamage_result.csv" > ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}_name_added_pydamage_result.csv

# no need to do it.
# mv ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}/pydamage_result.csv ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}_pydamage_result.csv

module unload bwa
module unload samtools
module unload bamtools
module unload pydamage

else
echo "Skipping DAMAGE PATTERN ANALYSIS."
fi

# KRAKEN FOR ASSEMBLIES

# NT database or custom database
DB="set-pathway-to-database"

# please do not change the level, otherwise, you will get nothing.
CONFIDENCE="0"

if [ "${RUN_KRAKEN}" = "YES" ]; then
module load kraken2
  if [[ ! -f ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${ID}_conf${CONFIDENCE}_contig.kraken ]]
then
    srun kraken2 --confidence ${CONFIDENCE} --db ${DB} ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${ID}/final.contigs.fa --threads ${CPU} --output ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${ID}_conf${CONFIDENCE}_contig.kraken \\
    --report ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${ID}_conf${CONFIDENCE}_contig.kraken.report
else
    echo "${ID} already exists on your filesystem [check "ref" directory]"
fi
echo ""
module unload kraken2/2.1.3
else
    echo "Skipping KRAKEN ANALYSIS."
