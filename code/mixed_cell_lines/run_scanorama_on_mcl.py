import scxx.preprocessing as pp
import scxx.plotting as pl
import scanorama
import numpy as np
import os
import scanpy as sc

test_result_path="./results/mixed_cell_lines/real/scanorama/"
import pandas as pd
genelist = pd.read_csv("./data/mixed_cell_lines/genelist.txt",header=None,index_col=0).index

if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)

r1="./data/mixed_cell_lines/293t.h5ad"
r2="./data/mixed_cell_lines/jurkat.h5ad"
r4="./data/mixed_cell_lines/mixed.h5ad"

adata_b1 = pp.read_sc_data(r1, batch_name="293t")
adata_b1 = adata_b1[:,genelist]
adata_b2 = pp.read_sc_data(r2, batch_name="jurkat")
adata_b2 = adata_b2[:,genelist]
adata_b4 = pp.read_sc_data(r4, batch_name="mixed")
adata_b4 = adata_b4[:,genelist]


datasets=[adata_b1,adata_b2,adata_b4]

integrated, corrected = scanorama.correct_scanpy(datasets, dimred=16, return_dimred=True)
scanorama_X=integrated[0]
adata_corrected=corrected[0]
adata_corrected.obs=datasets[0].obs
for i in np.arange(1,len(integrated)):
	scanorama_X=np.concatenate([scanorama_X,integrated[i]])
	adata=corrected[i]
	adata.obs=datasets[i].obs
	adata_corrected=adata_corrected.concatenate(adata,index_unique="-")

adata_corrected.raw=adata_corrected.copy()
adata_corrected.X=adata_corrected.X.todense()
adata_corrected.obsm["X_scanorama"]=scanorama_X
adata_corrected.obs.to_csv("./data/mixed_cell_lines/annotation_scanorama.txt",sep="\t")
adata_corrected.write(test_result_path+"output.h5ad")

