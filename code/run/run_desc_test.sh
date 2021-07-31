#!/bin/bash
source /home/jialu/packages/anaconda3/bin/activate descEnv
# python ./code/pancreas_scripts/run_desc_on_pancreas.py
# python ./code/bipolar_scripts/run_desc_on_bipolar.py
# python ./code/ctr_stim_pbmc_scripts/run_desc_on_ctr_stim_pbmc.py
python ./code/pbmcsca_scripts/run_desc_on_pbmcsca.py
python ./code/snrna_scripts/run_desc_on_mouse_brain.py