#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p debug
#SBATCH --mem 20g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/scmb_pluripotency_%j.out
#SBATCH -e ./jobs/scmb_pluripotency_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1 3 5;do
	for b1 in 64 128;do
		for b2 in 8 10 12 14 16 18 20;do
			python ./code/pluripotency_scripts/run_scxx_cvae_on_pluripotency.py CVAE $l ${b1} ${b2}
		done
	done
done
#Rscript --vanilla ./code/pluripotency_scripts/run_scxx_on_pluripotency.R
