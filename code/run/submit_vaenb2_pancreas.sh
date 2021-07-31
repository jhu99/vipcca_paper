#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p debug
#SBATCH --mem 30g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/vaenb2_pancreas_%j.out
#SBATCH -e ./jobs/vaenb2_pancreas_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 5;do
	for b1 in 32 64;do
			python ./code/pancreas_scripts/run_VAENB2_on_pancreas.py VAENB2 $l ${b1}
	done
done
