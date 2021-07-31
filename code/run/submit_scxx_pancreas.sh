#!/bin/bash
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -p control
#SBATCH --mem 20g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/scmb_pancreas_%j.out
#SBATCH -e ./jobs/scmb_pancreas_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 2;do
	for b1 in 64;do
		for b2 in 16 18;do
			python ./code/pancreas_scripts/run_scxx_on_pancreas.py CVAE $l ${b1} ${b2}
			#Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scxx_cvae/noscale_cvae_${l}_${b1}_${b2}/ $l_${b1}_${b2}
		done
	done
done
#Rscript --vanilla ./code/pancreas_scripts/run_evaluation_pancreas.R
