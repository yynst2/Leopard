#!/bin/bash

######## MOODS.sh ########
#This script is a wrapper around MOODS and will run MOODS using the motif, sequence, and cutoff specified

#INPUT: This script takes as input:

#INPUT_FASTA_DIR: directory of fastas
#INPUT_MOTIF_DIR: directory with motifs
#OUTPUT_DIR: directory were to write MOODS output
#PVAL: Pval cutoff to use
#MOODS_PARSER: Path to the parse_MOODS.py script
#TF_META: A tf meta file that associates a TF gene symbol to a motif
#TF: A list of TFs that are specific for the cell type of interest

#OUTPUT: 

#MOODS FILE: The output file from MOODS
#Predictions: BED file with predicted TF binding events
##########################

# Parameters
INPUT_FASTA=${1}
INPUT_MOTIF_DIR=${2}
OUTPUT_DIR=${3}
PVAL=${4}
TF_META=${5}
TF=${6}
PROJECT=${7}
INPUT_BINS=${8}
INPUT_BIN_TAG=${9}
METHOD=${10}

# Name Set up
INPUT_BASENAME=`basename ${1} .fa`
PROJECT_DIR=${OUTPUT_DIR}/${PROJECT}
MOODS_OUTPUT_FILENAME=${PROJECT_DIR}/MOODS/${PROJECT}_${TF}.mood
PROJECT_TF_META_FILE=${PROJECT_DIR}/${PROJECT}_${TF}_meta_file.tsv
TF_BEDS_DIR=${PROJECT_DIR}/BED
BINNED_PREDICTIONS_DIR=${PROJECT_DIR}/binned_predictions
MOODS_PREDICTIONS_DIR=${PROJECT_DIR}/MOODS
MOTIF_DIR=${PROJECT_DIR}/motifs/${TF}
OUT_TF_BED_FILENAME=${TF_BEDS_DIR}/${PROJECT}_${TF}.bed
OUT_GB_BED_FILENAME=${BINNED_PREDICTIONS_DIR}/${PROJECT}_${TF}_${INPUT_BIN_TAG}_${METHOD}.bed

echo "~~~~~~~~~~~~~ Run MOODS2BINS.sh ~~~~~~~~~~~~~"
echo "Input Fasta: " ${INPUT_FASTA}
echo "Motif Directory: " ${INPUT_MOTIF_DIR}
echo "Output Directory: " ${OUTPUT_DIR}
echo "TF BEDs Directory: " ${TF_BEDS_DIR}
echo "Binned Predictions Directory: " ${BINNED_PREDICTIONS_DIR}
echo "MOODS Predictions Directory: " ${MOODS_PREDICTIONS_DIR}
echo "MOODS Motifs used Directory: " ${MOTIF_DIR}
echo "MOODS Output Filename: " ${MOODS_OUTPUT_FILENAME}
echo "MOODS pvalue cutoff: " ${PVAL}
echo "Basename: " ${INPUT_BASENAME}
echo "TF Meta File: " ${TF_META}
echo "List of TFs to scan for: " ${TF}
echo "Project Name: " ${PROJECT}
echo "Project Directory: " ${PROJECT_DIR}
echo "Input Bin file: " ${INPUT_BINS}
echo "Input Bin tag: " ${INPUT_BIN_TAG}
echo "Method to perform on groupby of scores: " ${METHOD}
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~"

# Make the output directory if it does not exist
mkdir -p ${PROJECT_DIR}
mkdir -p ${MOTIF_DIR}
mkdir -p ${MOODS_PREDICTIONS_DIR}
mkdir -p ${TF_BEDS_DIR}
mkdir -p ${BINNED_PREDICTIONS_DIR}

echo "Motif_ID\tTF_Name" > ${PROJECT_TF_META_FILE}

# Get the meta information for the TFs in the list
grep -w ${TF} ${TF_META} | cut -f4,7 | sort | uniq >> ${PROJECT_TF_META_FILE}

# copy the motifs of interest into the output directory 
for i in $(tail -n +2 ${PROJECT_TF_META_FILE} | cut -f1 | sort | uniq);
    do
        cp ${INPUT_MOTIF_DIR}/${i}* ${MOTIF_DIR}/${i}
    done

# When using this script you must change into the directory with the motifs to prevent "an argument list too long" error
echo "Change to the motif directory"

cd ${MOTIF_DIR}

# Run MOODS
echo "Run MOODS"
moods-dna.py -m * -s ${INPUT_FASTA} -p ${PVAL} --batch -o ${MOODS_OUTPUT_FILENAME}

echo "Parse MOODS output into BED file"
# Once in a bed file format, bin the predictions and then groupby using bedtools. 
cut -d"," -f1,3,5,6 ${MOODS_OUTPUT_FILENAME} | \
tr ':' '\t' | tr '-' '\t' | tr ',' '\t' | \
awk '{ print $1"\t"$2+$4"\t"$2+$4+length($6)"\t"$5"\t"$1":"$2"-"$3}' | \
sort -k1,1 -k2,2n  > ${OUT_TF_BED_FILENAME}

echo "Bin the predictions and get the sum of predictions per bin"
bedtools intersect -wa -wb -sorted -nobuf -a ${INPUT_BINS} -b ${OUT_TF_BED_FILENAME} | \
bedtools groupby -g 1,2,3 -c 7 -o ${METHOD} > ${OUT_GB_BED_FILENAME}

echo "DONE!"