#!/bin/bash
## SBATCH -N 1
## SBATCH -n 40
## SBATCH -p debug
## SBATCH --mem 40g
## SBATCH -t 7-00:00
## SBATCH -o ./jobs/atacseq_pbmc_%j.out
## SBATCH -e ./jobs/atacseq_pbmc_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
# Run CVAE
for l in 3 4 5 10;do
	for b1 in 64 128;do
		for b2 in 10 16 32 64;do
			python ./code/atacseq_scripts/run_cvae_on_pbmc.py CVAE $l ${b1} ${b2}
		done
	done
done
# Rscript --vanilla ./code/pancreas_scripts/run_mnn_pancreas.R
#Rscript --vanilla ./code/pancreas_scripts/run_scalign_pancreas.R
# for p in ./results/pancreas/scxx_cvae_1e4/cvae_var2_30_*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae pancreas
# done
# for p in ./results/pancreas/scxx_cvae2/cvae2_30_*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae pancreas
# done
#kbet test, ARI for real data
#Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scvi/ scvi pancreas
#Rscript --vanilla ./code/base/kbet.R ./results/pancreas/desc/ desc pancreas
#Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scanorama/ scanorama pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/seurat/ seurat pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scalign/ scalign pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/liger/ liger pancreas
# Rscript --vanilla ./code/base/kbet.R ./results/pancreas/mnn/ mnn pancreas

