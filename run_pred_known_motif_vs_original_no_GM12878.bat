#PBS -l nodes=1:ppn=2:gpus=1
#PBS -l walltime=125:00:00
#PBS -l mem=80gb
#PBS -N pred_known_motif
#PBS -j oe
#PBS -m abe
#PBS -A PCCH0011

cd /users/PCCH0011/cch0017/PROJECTS/maxATAC/Leopard

source ~/.bashrc
source activate gpu


# known motif prediction
python Leopard_known_motif.py  -tf REST -te GM12878 -chr chr1 
python Leopard_known_motif.py  -tf CTCF -te GM12878 -chr chr1 
python Leopard_known_motif.py  -tf JUND -te GM12878 -chr chr1 

# original no GM12878 prediction
python Leopard_no_GM12878_seed=1.py  -tf REST -te GM12878 -chr chr1
python Leopard_no_GM12878_seed=1.py  -tf CTCF -te GM12878 -chr chr1
python Leopard_no_GM12878_seed=1.py  -tf JUND -te GM12878 -chr chr1

