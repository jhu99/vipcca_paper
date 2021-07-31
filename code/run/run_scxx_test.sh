#!/bin/bash
source /home/jialu/packages/anaconda3/bin/activate scEnv
# pancreas
# rm -rf ./results/pancreas/scxx_cvae/
# rm -rf ./results/pancreas/scxx_cvae2/
for l in 10;do
	for b1 in 64 128;do
		for b2 in 16 64 128;do
			python ./code/pancreas_scripts/run_scxx_on_pancreas.py CVAE $l ${b1} ${b2}
			#python ./code/pbmcsca_scripts/run_scxx_on_pbmcsca.py CVAE $l $b
			# python ./code/ctr_stim_pbmc_scripts/run_scxx_on_ctr_stim_pbmc.py CVAE $l ${b1} ${b2}
			# python ./code/ctr_stim_pbmc_scripts/run_scxx_on_ctr_stim_pbmc.py CVAE2 $l ${b1} ${b2}
		done
	done
done
#python ./code/pancreas_scripts/run_scxx_on_pancreas.py CVAE2
Rscript --vanilla ./code/pancreas_scripts/run_evaluation_pancreas.R
# # # pmbcsca
# rm -rf ./results/pbmcsca/scxx_cvae/
# rm -rf ./results/pbmcsca/scxx_cvae2/
# python ./code/pbmcsca_scripts/run_scxx_on_pbmcsca.py CVAE
# python ./code/pbmcsca_scripts/run_scxx_on_pbmcsca.py CVAE2
#Rscript --vanilla ./code/pbmcsca_scripts/run_evaluation_pbmcsca.R &
# # # bipolar
# rm -rf ./results/bipolar/scxx_cvae/
# rm -rf ./results/bipolar/scxx_cvae2/
# python ./code/bipolar_scripts/run_scxx_cvae_on_bipolar.py CVAE
# python ./code/bipolar_scripts/run_scxx_cvae_on_bipolar.py CVAE2
# # Rscript --vanilla ./code/bipolar_scripts/run_evaluation_bipolar.R &
# # # # kang's ctr_stim
# rm -rf ./results/ctr_stim_pbmc/scxx_cvae/
# rm -rf ./results/ctr_stim_pbmc/scxx_cvae2/
# python ./code/ctr_stim_pbmc_scripts/run_scxx_on_ctr_stim_pbmc.py CVAE
# python ./code/ctr_stim_pbmc_scripts/run_scxx_on_ctr_stim_pbmc.py CVAE2
# # Rscript --vanilla ./code/ctr_stim_pbmc_scripts/run_evaluation_ctr_stim.R  &
# # # # snrna
# rm -rf ./results/snrna/scxx_cvae/
# rm -rf ./results/snrna/scxx_cvae2/
# python ./code/snrna_scripts/run_scxx_on_mouse_brain.py CVAE
# python ./code/snrna_scripts/run_scxx_on_mouse_brain.py CVAE2
# Rscript --vanilla ./code/snrna_scripts/run_evaluation_on_mouse_brain.R  &
# querydata
# python ./code/querydata_scripts/run_scxx_on_querydata.py CVAE
# python ./code/querydata_scripts/run_scxx_on_querydata.py CVAE2
# hspc
# python ./code/hspc_scripts/run_scxx_on_hspc.py CVAE
# python ./code/hspc_scripts/run_scxx_on_hspc.py CVAE2
# python ./code/hspc_scripts/run_scxx_on_hspc_sub.py CVAE
# python ./code/hspc_scripts/run_scxx_on_hspc_sub.py CVAE2
