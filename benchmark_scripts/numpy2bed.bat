#BSUB -n 64
#BSUB -W 06:00
#BSUB -M 200000
#BSUB -o /data/miraldiNB/Tareian/std/NPY2BED_leopard.out
#BSUB -e /data/miraldiNB/Tareian/std/NPY2BED_leopard.err 

########## Body ##########
module load python3
module load gzip
module load bedtools 

N=8

(
for file in /data/miraldiNB/maxATAC/Leopard_output/SRX2717911_50bp_chr1/*.npy;
    do 
        ((i=i%N)); ((i++==0)) && wait

        sh /data/miraldiNB/maxATAC/bin/bin_Leopard_predictions.sh ${file} /data/miraldiNB/maxATAC/test_data/Leopard_output/SRX2717911_50bp_chr1/BED /data/miraldiNB/maxATAC/genome_inf/hg19/binned_genome/hg19_chr1_w200s200_blacklisted.bed w200 &
    done
)