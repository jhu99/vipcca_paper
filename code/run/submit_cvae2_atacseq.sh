#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p debug
#SBATCH --mem 20g
#SBATCH -t 4-00:00
#SBATCH -o ./jobs/cvae2_atacseq_%j.out
#SBATCH -e ./jobs/cvae2_atacseq_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 15 20;do
	for b1 in 10 12 14 16;do
		python ./code/atacseq_scripts/run_cvae2_on_atacseq.py CVAE2 $l ${b1}
	done
done
#Rscript --vanilla ./code/pancreas_scripts/run_evaluation_pancreas.R
