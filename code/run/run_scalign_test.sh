#!/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH -p debug
#SBATCH --cpus-per-task=9
#SBATCH -a 4
#SBATCH --mem 90g # memory pool for all cores
#SBATCH -t 5-10:00 # time (D-HH:MM)
#SBATCH -o ./jobs/scalign_%A_%a.out
#SBATCH -e ./jobs/scalign_%A_%a.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv_cpu
file=`ls ./code/*/*scalign* | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
echo $file
date +"%Y-%m-%d_%H:%M:%S"
Rscript --vanilla $file
date +"%Y-%m-%d_%H:%M:%S"
