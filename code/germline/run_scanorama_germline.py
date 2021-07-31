import vipcca.preprocessing as pp
import vipcca.plotting as pl
import scanorama
import numpy as np
import os
import scanpy as sc
import sys

part=str(sys.argv[1])
test_result_path="./results/germline/"+part+"/scanorama/"
rfile="./data/germline/"+part+"/"+part+".h5ad"
if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)

adata=pp.read_sc_data(rfile)
datasets=pp.split_object(adata,by="batch")

integrated, corrected = scanorama.correct_scanpy(datasets, dimred=16, return_dimred=True)
scanorama_X=integrated[0]
adata_corrected=corrected[0]
adata_corrected.obs=datasets[0].obs
for i in np.arange(1,len(integrated)):
	scanorama_X=np.concatenate([scanorama_X,integrated[i]])
	adata=corrected[i]
	adata.obs=datasets[i].obs
	adata_corrected=adata_corrected.concatenate(adata,index_unique=None)

adata_corrected.obs['batch']=adata_corrected.obs['batch'].cat.rename_categories({'0':"Guo",'1':"Li"})
adata_corrected.raw=adata_corrected.copy()
adata_corrected.X=adata_corrected.X.todense()
adata_corrected.obsm["X_scanorama"]=scanorama_X
adata_corrected.write(test_result_path+"output.h5ad")

