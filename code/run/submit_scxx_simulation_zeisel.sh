#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p debug
#SBATCH --mem 50g
#SBATCH -t 4-00:00
#SBATCH -o ./jobs/scmb_sim_zeisel_%j.out
#SBATCH -e ./jobs/scmb_sim_zeisel_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1;do
	for b1 in 64;do
		for b2 in 10;do
			python ./code/simulation/run_scxx_VAENB_simulation_zeisel.py VAENB $l ${b1} ${b2}
		done
	done
done

