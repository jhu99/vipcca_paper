import scxx.preprocessing as pp
import scxx.plotting as pl
import scanorama
import numpy as np
import os
import scanpy as sc

test_result_path="./results/mixed_cell_lines/permuted/scanorama/"
if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)

r1="./data/mixed_cell_lines/293t.h5ad"
r2="./data/mixed_cell_lines/jurkat.h5ad"
r3="./data/mixed_cell_lines/mixed_permuted.h5ad"
r4="./data/mixed_cell_lines/mixed.h5ad"

adata_b1 = pp.read_sc_data(r1, batch_name="293t", )
adata_b2 = pp.read_sc_data(r2, batch_name="jurkat")
adata_b3 = pp.read_sc_data(r3, batch_name="mixed_permutated")
#adata_b4 = pp.read_sc_data(r4, batch_name="mixed")


datasets=[adata_b1,adata_b2,adata_b3]

integrated, corrected = scanorama.correct_scanpy(datasets, dimred=16, return_dimred=True)
scanorama_X=integrated[0]
adata_corrected=corrected[0]
adata_corrected.obs=datasets[0].obs
for i in np.arange(1,len(integrated)):
	scanorama_X=np.concatenate([scanorama_X,integrated[i]])
	adata=corrected[i]
	adata.obs=datasets[i].obs
	adata_corrected=adata_corrected.concatenate(adata,index_unique=None)
	
import scanpy as scy
scy.pp.highly_variable_genes(adata_corrected,subset=True,n_top_genes=2000)
adata_corrected.raw=adata_corrected.copy()
scy.pp.scale(adata_corrected)

adata_corrected.obsm["X_scanorama"]=scanorama_X
adata_corrected.obs_names_make_unique()
adata_corrected.obs.to_csv("./data/mixed_cell_lines/annotation_scanorama_permuted.txt",sep="\t")
adata_corrected.write(test_result_path+"output.h5ad")

