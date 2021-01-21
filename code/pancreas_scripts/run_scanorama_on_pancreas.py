import scxx.preprocessing as pp
import scxx.plotting as pl
import scanorama
import os
import numpy as np
import scanpy as sc
from anndata import AnnData

test_input_file="./data/panc8/panc8.h5ad"
test_result_path="./results/pancreas/scanorama/"
import pandas as pd
genelist = pd.read_csv("./data/panc8/genelist.txt",header=None,index_col=0).index
if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)

adata = pp.read_sc_data(test_input_file)
adata = adata[:,genelist]
datasets=pp.split_object(adata,by="dataset")

integrated, corrected = scanorama.correct_scanpy(datasets, return_dimred=True, dimred=16)

scanorama_X=integrated[0]
adata_corrected=corrected[0]
adata_corrected.obs=datasets[0].obs
for i in np.arange(1,len(integrated)):
	scanorama_X=np.concatenate([scanorama_X,integrated[i]])
	adata_i=corrected[i]
	adata_i.obs=datasets[i].obs
	adata_corrected=adata_corrected.concatenate(adata_i,index_unique=None)
adata_corrected.raw=adata_corrected.copy()
adata_corrected.X=adata_corrected.X.todense()
adata_corrected.obsm["X_scanorama"]=scanorama_X
adata_corrected.write(test_result_path+"corrected.h5ad")
# 
