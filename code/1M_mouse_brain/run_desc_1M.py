import scanpy as sc
#import scxx.preprocessing as pp
import desc
import os
import pandas as pd

test_result_path="./results/1M_mouse_brain/desc/"

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

genelist = pd.read_csv("./data/scanorama_data/data/mouse_brain/genelist_vipcca.txt",header=None,index_col=0).index

for i in range(len(data_names)):
        name=data_names[i]
        ann = sc.read_h5ad("./data/scanorama_data/"+name+"/data.h5ad")
        ann.obs[".batch"]=i
        if i==0:
        	adata=ann
        else:
        	adata=adata.concatenate(ann)
        
# 
sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e6)
sc.pp.log1p(adata)
adata=adata[:,genelist]
adataw=adata.copy()
sc.pp.scale(adata,zero_center=True, max_value=3)

adata=desc.train(adata,
        dims=[adata.shape[1],64,32,16],
        tol=0.005,
        n_neighbors=10,
        batch_size=256,
        louvain_resolution=0.3,# not necessarily a list, you can only set one value, like, louvain_resolution=1.0
        save_dir=str(test_result_path),
        do_tsne=False,
        learning_rate=200, # the parameter of tsne
        use_GPU=True,
        num_Cores=24, #for reproducible, only use 1 cpu
        num_Cores_tsne=10,
        save_encoder_weights=True,
        save_encoder_step=3,# save_encoder_weights is False, this parameter is not used
        use_ae_weights=False,
        do_umap=False)
adata.raw=adataw
adata.write(test_result_path+"output.h5ad")

#adata=sc.read_h5ad(test_result_path+"output.h5ad")
import pandas as pd
obs=pd.read_csv("./results/1M_mouse_brain/barcode_all.txt",header=None)
obs_selected=pd.read_csv("./results/1M_mouse_brain/barcode_selected.txt",header=None)
adata.obs.index = obs.iloc[:,0]
adata.obs.index.rename("barcode",inplace=True)
adata_sub=adata[obs_selected.iloc[:,0],:]
adata_sub.write(test_result_path+"output_sub.h5ad")
