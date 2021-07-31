#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p control
#SBATCH --mem 50g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/scmb_bipolar_%j.out
#SBATCH -e ./jobs/scmb_bipolar_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1 3 5;do
	for b1 in 64 128;do
		for b2 in 8 12 16 18 20;do
			python ./code/bipolar_scripts/run_scxx_cvae_on_bipolar.py CVAE $l ${b1} ${b2}
		done
	done
done
#Rscript --vanilla ./code/bipolar_scripts/run_evaluation_bipolar.R
