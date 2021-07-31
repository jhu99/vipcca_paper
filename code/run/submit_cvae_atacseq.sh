#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p control
#SBATCH --mem 20g
#SBATCH -t 4-00:00
#SBATCH -o ./jobs/cvae_atacseq_%j.out
#SBATCH -e ./jobs/cvae_atacseq_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 5 10;do
	for b1 in 64;do
		for b2 in 10 20;do
			python ./code/atacseq_scripts/run_cvae_on_atacseq.py CVAE $l ${b1} ${b2}
		done
	done
done
#Rscript --vanilla ./code/pancreas_scripts/run_evaluation_pancreas.R
