import scxx.preprocessing as pp
import scxx.plotting as pl
import scanorama
import os
import numpy as np
import scanpy as sc
from anndata import AnnData

np.random.seed(0)

NAMESPACE = 'mouse_brain'
BATCH_SIZE = 1000
result_dir="./results/1M_mouse_brain/scanorama/"

data_names = [
    'data/mouse_brain/nuclei',
    'data/mouse_brain/dropviz/Cerebellum_ALT',
    'data/mouse_brain/dropviz/Cortex_noRep5_FRONTALonly',
    'data/mouse_brain/dropviz/Cortex_noRep5_POSTERIORonly',
    'data/mouse_brain/dropviz/EntoPeduncular',
    'data/mouse_brain/dropviz/GlobusPallidus',
    'data/mouse_brain/dropviz/Hippocampus',
    'data/mouse_brain/dropviz/Striatum',
    'data/mouse_brain/dropviz/SubstantiaNigra',
    'data/mouse_brain/dropviz/Thalamus',
]
import pandas as pd
genelist = pd.read_csv("./data/scanorama_data/data/mouse_brain/genelist_vipcca.txt",header=None,index_col=0).index


datasets=[]
for i in range(len(data_names)):
	name=data_names[i]
	ann = pp.read_sc_data("./data/scanorama_data/"+name+"/data.h5ad",batch_name=str(i))
	ann=ann[:,genelist]
	ann.write("./data/scanorama_data/"+name+"/data_subset.h5ad")
	datasets.append(ann)

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
adata_corrected.obs_names_make_unique()

# 1,094,150
adata_corrected.write(result_dir+"output.h5ad")
