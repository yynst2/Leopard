#!/bin/bash

###### bin_Leopard_prediction.sh #########

# Input
# 1 ) NP_array
# 2 ) Output directory to write results to
# 3 ) The genome bins bed file
# 4 ) string to use for naming the files based on bins

Leopard_NP_array=${1}
OUT_DIR=${2}
BINS=${3}
window_tag=${4}

# Make the file name for the 1 bp resolution Leopard results
Leopard_BED=`basename ${Leopard_NP_array} .npy`_1bp.bed

# This is the python script used to convert NP array to bed
Leopard2BED_script="/data/miraldiNB/maxATAC/bin/numpy2bed.py"

mkdir -p ${OUT_DIR}

# Change to target output directory
cd ${OUT_DIR}

# Parse NP to BED
python3 ${Leopard2BED_script} -i ${Leopard_NP_array}

# Intersect the predictions with the bins and then groupby the 7th column (In my example; change according to needs) and sum the results
bedtools intersect -wa -wb -nobuf -sorted -a ${BINS} -b ${Leopard_BED} | bedtools groupby -g 1,2,3 -c 7 -o sum > `basename ${Leopard_BED} _1bp.bed`_${window_tag}_sum.bed

