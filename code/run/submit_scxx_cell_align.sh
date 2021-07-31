#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p debug
#SBATCH --mem 20g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/scmb_cell_align_%j.out
#SBATCH -e ./jobs/scmb_cell_align_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1 3 5;do
	for b1 in 64 128;do
		for b2 in 8 10 12 14 16 18 20;do
			python ./code/cell_align/run_scxx_cvae_on_dc.py CVAE $l ${b1} ${b2}
		done
	done
done
#Rscript --vanilla ./code/pluripotency_scripts/run_scxx_on_pluripotency.R
