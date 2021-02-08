#!/bin/bash
source /home/${USER}/.bashrc
RESPATH="./results/1M_mouse_brain/seurat"
#source /home/jialu/packages/anaconda3/bin/activate scEnv
# Run CVAE
# cmd="python ./code/1M_mouse_brain/run_scxx_cvae_on_1M.py CVAE 2 64 16"
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >log.stdout 2>log.stderr & } && echo $!);

# Run seurat
# cmd="Rscript --vanilla ./code/1M_mouse_brain/run_seurat_1M.R"
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# # Run liger
# RESPATH="./results/1M_mouse_brain/liger"
# cmd="Rscript --vanilla ./code/1M_mouse_brain/run_liger_1M.R"
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# Run mnn
# RESPATH="./results/1M_mouse_brain/mnn"
# cmd="Rscript --vanilla ./code/1M_mouse_brain/run_mnn_1M.R"
# mkdir ${RESPATH}
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# Run scVI
# source /home/jialu/packages/anaconda3/bin/activate scviEnv
# RESPATH="./results/1M_mouse_brain/scvi"
# cmd="python ./code/1M_mouse_brain/run_scvi_1M.py"
# mkdir ${RESPATH}
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# # Run DESC
# source /home/jialu/packages/anaconda3/bin/activate descEnv
# RESPATH="./results/1M_mouse_brain/desc"
# cmd="python ./code/1M_mouse_brain/run_desc_1M.py"
# mkdir ${RESPATH}
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# Run scalign
# source /home/jialu/packages/anaconda3/bin/activate scEnv
# RESPATH="./results/1M_mouse_brain/scalign"
# cmd="Rscript --vanilla ./code/1M_mouse_brain/run_scalign_1M.R"
# mkdir ${RESPATH}
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} >${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# RESPATH="./results/1M_mouse_brain/scanorama"
# source /home/jialu/packages/anaconda3/bin/activate scanoramaEnv
# mkdir ./results/1M_mouse_brain/scanorama
# cmd="python ./code/1M_mouse_brain/run_scanorama_1M.py"
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} > ${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

# RESPATH="./results/1M_mouse_brain/scanorama_sketched"
# source /home/jialu/packages/anaconda3/bin/activate scanoramaEnv
# mkdir ./results/1M_mouse_brain/scanorama_sketched
# cmd="python ./code/1M_mouse_brain/mouse_brain_sketched.py"
# echo "# ${cmd}"
# WATCHED_PID=$({ ${cmd} > ${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

RESPATH="./results/1M_mouse_brain/harmony"
source /home/jialu/packages/anaconda3/bin/activate scEnv
mkdir ./results/1M_mouse_brain/harmony
cmd="Rscript --vanilla ./code/1M_mouse_brain/run_harmony_1M.R"
echo "# ${cmd}"
WATCHED_PID=$({ ${cmd} > ${RESPATH}/log.stdout 2>${RESPATH}/log.stderr & } && echo $!);

top -b -d 1 -p $WATCHED_PID | awk -v OFS="," '$1+0>0 {print $9,$10; fflush() }'