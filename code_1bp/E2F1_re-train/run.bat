#PBS -l nodes=1:ppn=2:gpus=1:default
#PBS -l walltime=125:00:00
#PBS -l mem=50gb
#PBS -N E2F1
#PBS -j oe
#PBS -m abe


cd /users/PCCH0011/cch0017/PROJECTS/maxATAC/Leopard/code_1bp/E2F1_re-train


source ~/.bashrc
source activate gpu

python train.py -tf E2F1 -tr  HeLa-S3  -vali GM12878   -par 1
