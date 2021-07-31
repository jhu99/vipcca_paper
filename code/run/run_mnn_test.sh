#!/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH -p control 
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 5-10:00 # time (D-HH:MM)
#SBATCH -o ./jobs/mnn_%A_%a.out
#SBATCH -e ./jobs/mnn_%A_%a.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
# file=`ls ./code/[^r]*/*mnn* | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
# echo $file
# date +"%Y-%m-%d_%H:%M:%S"
# Rscript --vanilla $file
# date +"%Y-%m-%d_%H:%M:%S"
Rscript --vanilla ./code/pbmcsca/run_mnn_pbmcsca.R 
Rscript --vanilla ./code/base/kbet.R ./results/pbmcsca/mnn/ mnn pbmcsca
