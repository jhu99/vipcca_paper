#!/bin/bash
# #SBATCH -N 1
# #SBATCH -n 15
# #SBATCH -p control
# #SBATCH --mem 50g
# #SBATCH -t 4-00:00
# #SBATCH -o ./jobs/scmb_sim_visp_%j.out
# #SBATCH -e ./jobs/scmb_sim_visp_%j.err
# source /home/${USER}/.bashrc
# source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1 3 5;do
	for b1 in 8 64;do
		for b2 in 8 64;do
			python ./code/simulation/run_scxx_cvae_simulation_VISp.py CVAE $l ${b1} ${b2}
		done
	done
done