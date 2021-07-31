#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p control
#SBATCH --mem 10g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/scmb_pbmcsca_%j.out
#SBATCH -e ./jobs/scmb_pbmcsca_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
for l in 1;do
	for b1 in 64 128;do
		for b2 in 10 14;do
			python ./code/pbmcsca_scripts/run_scxx_on_pbmcsca.py CVAE $l ${b1} ${b2}
			Rscript --vanilla ./code/base/kbet.R ./results/pbmcsca/scxx_cvae/scexp_cvae_${l}_${b1}_${b2}/ ${l}${b1}${b2}
		done
	done
done

