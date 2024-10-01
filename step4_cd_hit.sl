#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: CD-HIT gene cluster to construct non-redundant gene catalog
#===========================================================================


WORK=${PWD}
OUTDIR="output"
OUT_PRODIGAL="out.prodigal"
OUT_METAEUK="out.metaeuk"
OUT_TIARA="out.tiara"

# Example for prokGAP
IN="redundant_PROKGAP_protein.faa"
OUT="non_redundant_PROKGAP_protein.faa"

cd-hit -i ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${IN} -o ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${OUT} -n 5 -c 0.95 -G 0 -M 0 -d 0 -aS 0.85
