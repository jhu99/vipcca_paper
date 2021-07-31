#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p control
#SBATCH --mem 50g
#SBATCH -t 4-00:00
#SBATCH -o ./jobs/vaenb2_sim_visp_%j.out
#SBATCH -e ./jobs/vaenb2_sim_visp_%j.err
# source /home/${USER}/.bashrc
# source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1;do
	for b1 in 32 64;do
		python ./code/simulation/run_scxx_VAENB2_simulation_VISp.py VAENB2 $l ${b1}
	done
done

