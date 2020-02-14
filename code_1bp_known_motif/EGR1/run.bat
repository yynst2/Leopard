#PBS -l nodes=1:ppn=2:gpus=1
#PBS -l walltime=125:00:00
#PBS -l mem=40gb
#PBS -N E2F1
#PBS -j oe
#PBS -m abe


cd /users/PCCH0011/cch0017/PROJECTS/maxATAC/Leopard/code_1bp_known_motif/E2F1


source ~/.bashrc
source activate gpu

python /users/PCCH0011/cch0017/PROJECTS/maxATAC/Leopard/unet_known_motif_structure/train.py \
-tf E2F1 \
-tr  HeLa-S3  \
-vali GM12878   \
-par 1
