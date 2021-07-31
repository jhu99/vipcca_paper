#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -p control
#SBATCH --mem 30g
#SBATCH -t 4-00:00
#SBATCH -o ./jobs/cvae2_pbmcsca_%j.out
#SBATCH -e ./jobs/cvae2_pbmcsca_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1 3 5;do
	for b1 in 128;do
			#python ./code/pbmcsca_scripts/run_cvae2_on_pbmcsca.py CVAE2 $l ${b1}
			Rscript --vanilla ./code/base/kbet.R ./results/pbmcsca/scxx_cvae2/scexp_cvae2_${l}_${b1}/ ${l}${b1}
	done
done
