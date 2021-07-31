#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p control
#SBATCH --mem 30g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/cvae2_pancreas_%j.out
#SBATCH -e ./jobs/cvae2_pancreas_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 3 4;do
	for b1 in 8 16 64 128;do
			python ./code/pancreas_scripts/run_cvae2_on_pancreas.py CVAE2 $l ${b1}
			Rscript --vanilla ./code/base/kbet.R ./results/pancreas/scxx_cvae2/noscale_var1_${l}_${b1}/ ${l}_${b1}
	done
done
