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
#TF_LIST: A list of TFs that are specific for the cell type of interest

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
TF_LIST=${6}
PROJECT=${7}

INPUT_BASENAME=`basename ${1} .fa`
MOODS_PARSER="/Users/caz3so/workspaces/miraldilab/team/Tareian/python/MOODS2TFBEDS_parallel.py"
PROJECT_DIR=${OUTPUT_DIR}/${PROJECT}

echo "Input Fasta: " ${INPUT_FASTA}
echo "Motif Directory: " ${INPUT_MOTIF_DIR}
echo "Output Directory: " ${OUTPUT_DIR}
echo "MOODS pvalue cutoff: " ${PVAL}
echo "Basename: " ${INPUT_BASENAME}
echo "Python script to parse MOODS output: " ${MOODS_PARSER}
echo "TF Meta File: " ${TF_META}
echo "List of TFs to scan for: " ${TF_LIST}
echo "Project Name: " ${PROJECT}
echo "Project Directory: " ${PROJECT_DIR}

# Make the output directory if it does not exist
mkdir -p ${PROJECT_DIR}
mkdir -p ${PROJECT_DIR}/motifs
mkdir -p ${PROJECT_DIR}/MOODS
mkdir -p ${PROJECT_DIR}/TF_BEDS

echo "Motif_ID\tTF_Name" > ${PROJECT_DIR}/${PROJECT}_meta_file.tsv

# Get the meta information for the TFs in the list
grep -w -f ${TF_LIST} ${TF_META} | cut -f4,7 | sort | uniq >> ${PROJECT_DIR}/${PROJECT}_meta_file.tsv

# copy the motifs of interest into the output directory 
for i in $(tail -n +2 ${PROJECT_DIR}/${PROJECT}_meta_file.tsv | cut -f1 | sort | uniq);
    do
        cp ${INPUT_MOTIF_DIR}/${i}* ${PROJECT_DIR}/motifs/${i}
    done

# When using this script you must change into the directory with the motifs to prevent "an argument list too long" error
echo "Change to the motif directory"

cd ${PROJECT_DIR}/motifs

# Run MOODS
echo "Run MOODS"
moods-dna.py -m * -s ${INPUT_FASTA} -p ${PVAL} --batch --sep "|" -o ${PROJECT_DIR}/MOODS/${PROJECT}.mood

# The MOODS output file is not in the format that we want it in so we have to parse it. 
echo "Parse the MOODS output and write a bed file for each TF"
python ${MOODS_PARSER} -i ${PROJECT_DIR}/MOODS/${PROJECT}.mood -T ${PROJECT_DIR}/${PROJECT}_meta_file.tsv -K ${TF_LIST} -o ${PROJECT_DIR}/TF_BEDS

echo "DONE!"
