#!/bin/bash
#
source activate vipcca
for l in 10 20 30 40 50;do
  for gender in {"male","female"}; do
    for rp in {1..5};do
      python ./code/germline/run_scxx_cvae_on_germline.py CVAE ${l} 64 16 ${gender} ${rp}
		done
	done
done
#
# source activate desc
# python ./code/germline/run_desc_on_germline.py female
# python ./code/germline/run_desc_on_germline.py male
# source activate scanorama
# python ./code/germline/run_scanorama_germline.py female
# python ./code/germline/run_scanorama_germline.py male
# Rscript ./code/germline/run_mnn_germline.R female
# Rscript ./code/germline/run_mnn_germline.R male
# Rscript ./code/germline/run_scalign_on_germline.R female
# Rscript ./code/germline/run_scalign_on_germline.R male
# Rscript ./code/germline/run_seurat_germline.R female
# Rscript ./code/germline/run_seurat_germline.R male
# Rscript ./code/germline/run_harmony_germline.R female
# Rscript ./code/germline/run_harmony_germline.R male
# # Rscript ./code/germline/run_liger_germline.R female
# # Rscript ./code/germline/run_liger_germline.R male
# Rscript ./code/germline/run_uncorrected_germline.R female
# Rscript ./code/germline/run_uncorrected_germline.R male



