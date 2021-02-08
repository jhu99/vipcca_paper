# !/bin/bash
# # SBATCH -N 1
# # SBATCH -n 5
# # SBATCH -p compute
# # SBATCH --mem 20g
# # SBATCH -t 2-00:00
# # SBATCH -o ./jobs/mcl_%j.out
# # SBATCH -e ./jobs/mcl_%j.err
source /home/${USER}/.bashrc
source /home/jialu/packages/anaconda3/bin/activate scEnv
## Run CVAE
# for l in 2;do
# 	for b1 in 128;do
# 		for b2 in 8 10 12 14 16;do
# 			python ./code/mixed_cell_lines/run_scxx_cvae_on_mcl_permuted.py CVAE $l ${b1} ${b2}
# 			# python ./code/mixed_cell_lines/run_scxx_cvae_on_mcl.py CVAE $l ${b1} ${b2}
# 		done
# 	done
# done
# for p in ./results/mixed_cell_lines/real/vipcca/cvae*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae mcl
# done
# for p in ./results/mixed_cell_lines/permuted/vipcca/cvae*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae mcl
# done

#Rscript --vanilla ./code/mixed_cell_lines/run_seurat_mcl.R
#Rscript --vanilla ./code/mixed_cell_lines/run_scalign_on_mcl.R
#Rscript --vanilla ./code/mixed_cell_lines/run_mnn_mcl.R
#Rscript --vanilla ./code/mixed_cell_lines/run_liger_mcl.R
#python ./code/mixed_cell_lines/run_scxx_cvae_on_mcl.py CVAE 2 64 10
#source /home/jialu/packages/anaconda3/bin/activate scviEnv
#python ./code/mixed_cell_lines/run_scvi_mcl.py
#python ./code/mixed_cell_lines/run_desc_on_mcl.py
### Run scanorama
#source /home/jialu/packages/anaconda3/bin/activate scanoramaEnv
#python ./code/mixed_cell_lines/run_scanorama_on_mcl.py
# python ./code/mixed_cell_lines/run_scanorama_on_mcl_permuted.py

#kbet test, ARI for real data
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/vipcca/cvae_2_64_14/ cvae mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/scanorama/ scanorama mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/desc/ desc mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/scalign/ scalign mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/mnn/ mnn mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/liger/ liger mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/seurat/ seurat mcl
Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/real/harmony/ harmony mcl
#Rscript --vanilla ./code/mixed_cell_lines/run_seurat_mcl.R
#Rscript --vanilla ./code/mixed_cell_lines/run_scalign_on_mcl.R
#Rscript --vanilla ./code/mixed_cell_lines/run_mnn_mcl.R
#Rscript --vanilla ./code/mixed_cell_lines/run_liger_mcl.R
#python ./code/mixed_cell_lines/run_scxx_cvae_on_mcl.py CVAE 2 64 10
#source /home/jialu/packages/anaconda3/bin/activate scviEnv
#python ./code/mixed_cell_lines/run_scvi_mcl.py
#python ./code/mixed_cell_lines/run_desc_on_mcl.py
### Run scanorama
#source /home/jialu/packages/anaconda3/bin/activate scanoramaEnv
#python ./code/mixed_cell_lines/run_scanorama_on_mcl.py
# python ./code/mixed_cell_lines/run_scanorama_on_mcl_permuted.py

# for p in ./results/mixed_cell_lines/permuted/vipcca/cvae*; do
# 	echo "$p/"
# 	Rscript --vanilla ./code/base/kbet.R ${p}/ cvae mcl
# done
# #kbet test, ARI for permuted data
# Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/permuted/scxx_cvae/cvae_2_64_8/ cvae mcl
# Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/permuted/liger/ liger mcl
# Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/permuted/seurat/ seurat mcl
# Rscript --vanilla ./code/base/kbet.R ./results/mixed_cell_lines/permuted/scanorama/ scanorama mcl
# #Plots