#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p compute
#SBATCH --mem 80g
#SBATCH -t 2-00:00
#SBATCH -o ./jobs/scanorama_pbmcsca_%j.out
#SBATCH -e ./jobs/scanorama_pbmcsca_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
# file=`ls ./code/[^r]*/*scanorama*`
# for scr in $file
# do
# python $scr &
# done
# python ./code/ctr_stim_pbmc_scripts/run_scanorama_on_ctr_stim_pbmc.py
# python ./code/pbmcsca_scripts/run_scanorama_on_pbmcsca.py
# python ./code/snrna_scripts/run_scanorama_on_mouse_brain.py
Rscript --vanilla ./code/base/kbet.R ./results/pbmcsca/scanorama/ scanorama pbmscsca