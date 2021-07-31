#!/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH --cpus-per-task=4
#SBATCH -p control
#SBATCH -a 0-1
#SBATCH --mem 50g # memory pool for all cores
#SBATCH -t 1-10:00 # time (D-HH:MM)
#SBATCH -o ./jobs/evaluation_%A_%a.out
#SBATCH -e ./jobs/evaluation_%A_%a.err

source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv_cpu
declare -a eval_script=(
                 "./code/pancreas_scripts/run_evaluation_pancreas.R"
                 "./code/pbmcsca_scripts/run_evaluation_pbmcsca.R"
                 "./code/bipolar_scripts/run_evaluation_bipolar.R" 
                 "./code/ctr_stim_pbmc_scripts/run_evaluation_ctr_stim.R"
                 "./code/snrna_scripts/run_evaluation_on_mouse_brain.R"
                )
declare -a scxx_script=(
                 "./code/pancreas_scripts/run_scxx_on_pancreas.py"
                 "./code/pbmcsca_scripts/run_scxx_on_pbmcsca.py"
                 "./code/bipolar_scripts/run_scxx_cvae_on_bipolar.py"
                 "./code/ctr_stim_pbmc_scripts/run_scxx_on_ctr_stim_pbmc.py"
                 "./code/snrna_scripts/run_scxx_on_mouse_brain.py"
                )
# python ${scxx_script[$SLURM_ARRAY_TASK_ID]} CVAE
# python ${scxx_script[$SLURM_ARRAY_TASK_ID]} CVAE2
Rscript --vanilla ${eval_script[$SLURM_ARRAY_TASK_ID]} 

###
# rm ./results/bipolar/scxx_cvae*/*.rds
# rm results/ctr_stim_pbmc/scxx_cvae*/*rds
# rm ./results/pancreas/scxx_cvae*/*rds
# rm ./results/pbmcsca/scxx_cvae*/*rds
# rm ./results/snrna/scxx_cvae*/*rds
