#!/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH --cpus-per-task=4
#SBATCH -p control
#SBATCH -a 3-4
#SBATCH --mem 60g # memory pool for all cores
#SBATCH -t 1-10:00 # time (D-HH:MM)
#SBATCH -o ./jobs/liger_%A_%a.out
#SBATCH -e ./jobs/liger_%A_%a.err

source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv_cpu
declare -a file=(
                 "./code/bipolar_scripts/run_liger_bipolar.R" 
                 "./code/pancreas_scripts/run_liger_pancreas.R" 
                 "./code/ctr_stim_pbmc_scripts/run_liger_ctr_stim.R"
                 "./code/pbmcsca_scripts/run_liger_pbmcsca.R"
                 "./code/snrna_scripts/run_liger_mouse_brain.R"
                )
Rscript --vanilla ${file[$SLURM_ARRAY_TASK_ID]} 

###
# rm ./results/bipolar/scxx_cvae*/*.rds
# rm results/ctr_stim_pbmc/scxx_cvae*/*rds
# rm ./results/pancreas/scxx_cvae*/*rds
# rm ./results/pbmcsca/scxx_cvae*/*rds
# rm ./results/snrna/scxx_cvae*/*rds
