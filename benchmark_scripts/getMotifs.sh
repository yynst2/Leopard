############ getMotifs.sh ############
# This script will take in a:
# INPUT 1: Directory of CISBP motifs
# INPUT 2: CISBP Meta File
# INPUT 3: Output directory for motifs
# INPUT 4: List of TFs you want to find motifs for

# This script will output a directory with motifs for every TF you want to analyze (TF List)

######################################

IN_MOTIF_DIR=${1}
META_FILE=${2}
OUTPUT_DIR=${3}
TF_LIST=${4}
PROJECT=${5}

# Make the output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

grep -w -f ${TF_LIST} ${META_FILE} | cut -f4,7 | sort | uniq >> ${PROJECT}_meta_file.txt

for i in $(grep -w -f ${TF_LIST} ${META_FILE} | cut -f4 | sort | uniq);
    do
        cp ${IN_MOTIF_DIR}/$i ${OUTPUT_DIR}/${i}
    done
