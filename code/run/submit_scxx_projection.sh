#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p debug
#SBATCH --mem 50g
#SBATCH -t 4-00:00
#SBATCH -o ./jobs/scmb_projection_%j.out
#SBATCH -e ./jobs/scmb_projection_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1 3 5;do
	for b1 in 8 12 16;do
		for b2 in 64 128;do
			python ./code/starmap_scripts/run_scxx_cvae_on_projection.py CVAE $l ${b1} ${b2}
		done
	done
done

