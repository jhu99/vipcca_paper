#!/bin/bash
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p debug
#SBATCH --mem 40g
#SBATCH -t 2-00:00
#SBATCH --array 0-7
#SBATCH -o ./jobs/pancreas_%A_%a.out
#SBATCH -e ./jobs/pancreas_%A_%a.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
# Run CVAE
# for l in 30;do
# 	for b1 in 10 16 32 64 128;do
# 		for b2 in 10 16 32 64 128;do
# 			python ./code/pancreas_scripts/run_scxx_on_pancreas.py CVAE $l ${b1} ${b2}
# 		done
# 		# python ./code/pancreas_scripts/run_cvae2_on_pancreas.py CVAE2 $l ${b1}
# 	done
# done
# Rscript --vanilla ./code/pancreas_scripts/run_mnn_pancreas.R
#Rscript --vanilla ./code/pancreas_scripts/run_scalign_pancreas.R
# for p in ./results/pancreas/scxx_cvae_1e4/cvae_var2_30_*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae pancreas
# done
# for p in ./results/pancreas/scxx_cvae/cvae_var2_$SLURM_ARRAY_TASK_ID_*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae pancreas
# done
#kbet test, ARI for real data
par=(
"./results/pancreas/scxx_cvae/cvae_var2_10_128_16/ cvae pancreas"
"./results/pancreas/desc/ desc pancreas"
"./results/pancreas/scanorama/ scanorama pancreas"
"./results/pancreas/seurat/ seurat pancreas"
"./results/pancreas/scalign/ scalign pancreas"
"./results/pancreas/liger/ liger pancreas"
"./results/pancreas/mnn/ mnn pancreas"
"./results/pancreas/harmony/ harmony pancreas")
echo ${par[$SLURM_ARRAY_TASK_ID]}
Rscript --vanilla ./code/base/kbet.R ${par[$SLURM_ARRAY_TASK_ID]}

# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scxx_cvae/cvae_var2_10_128_16/ cvae pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/desc/ desc pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scanorama/ scanorama pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/seurat/ seurat pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scalign/ scalign pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/liger/ liger pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/mnn/ mnn pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/harmony/ harmony pancreas

