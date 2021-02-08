import os
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from scvi.dataset import CortexDataset, RetinaDataset
from scvi.models import *
from scvi.inference import UnsupervisedTrainer
from scvi.inference.autotune import auto_tune_scvi_model
import scvi
from scvi.dataset import AnnDatasetFromAnnData
from scvi.models.vae import VAE
from typing import Tuple
import scxx.preprocessing as pp
import torch
from sklearn.metrics import mean_squared_error
from sklearn.metrics import median_absolute_error
import matplotlib
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy.sparse import csr_matrix
import seaborn as sns
matplotlib.use('Agg')
sns.set(style='white', rc={'figure.figsize':(5,5), 'figure.dpi':150})
fontsize=10
params = {'legend.fontsize': fontsize,
          'figure.figsize': (8, 8),
          'figure.dpi': 80,
         'axes.labelsize': fontsize,
         'axes.titlesize':fontsize,
         'xtick.labelsize':fontsize,
         'ytick.labelsize':fontsize}
plt.rcParams.update(params)

lr=1e-3
show_plot=False

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


for i in range(len(data_names)):
	name=data_names[i]
	ann = pp.read_sc_data("./data/scanorama_data/"+name+"/data.h5ad",batch_name=str(i))
	if i==0:
		annAll=ann
	else:
		annAll=annAll.concatenate(ann)
		
import pandas as pd
obs=pd.read_csv("./results/1M_mouse_brain/barcode_all.txt",header=None)
annAll.obs.index=obs.iloc[:,0]
annAll.obs.index.rename("barcode",inplace=True)
genelist = pd.read_csv("./data/scanorama_data/data/mouse_brain/genelist_vipcca.txt",header=None,index_col=0).index
annAll = annAll[:,genelist]
import scanpy as sc
sc.pp.filter_cells(annAll, min_counts=1)

test_result_path="./results/1M_mouse_brain/scVI/"
model_file=test_result_path + "model.pkl"
if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)
scviDataset = AnnDatasetFromAnnData(annAll, batch_label="batch")

use_cuda=True
use_batches=True
n_epochs = 500
mse=np.zeros(10)
mae=np.zeros(10)
k=0

vae = VAE(scviDataset.nb_genes, 
					n_batch=scviDataset.n_batches * use_batches,
					n_hidden=128, 
					n_latent=16, 
					n_layers=1,
					dispersion='gene')
trainer = UnsupervisedTrainer(vae,
                              scviDataset,
                              train_size=1.0,
                              use_cuda=use_cuda
                              )
if os.path.isfile(model_file):
	trainer.model.load_state_dict(torch.load(model_file))
	trainer.model.eval()
else:
	trainer.train(n_epochs=n_epochs, lr=lr)
	torch.save(trainer.model.state_dict(), model_file)

full = trainer.create_posterior(trainer.model, scviDataset, indices=np.arange(len(scviDataset)))
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

cell_freq=full.sequential().get_sample_scale()
adata_t=annAll.copy()
adata_t.X=cell_freq
cell_freq_copy=cell_freq.copy()
cell_freq_copy[cell_freq<1e-6]=0
adata_t.raw=AnnData(X=csr_matrix(cell_freq_copy))
adata_t.obsm['X_scvi']=latent
adata_t.write(test_result_path+"output.h5ad")

obs_selected=pd.read_csv("./results/1M_mouse_brain/barcode_selected_2.txt",header=None)
ind_selected=adata_t.obs.index.intersection(obs_selected.iloc[:,0])
res_sub=adata_t[ind_selected,:]
res_sub.write(test_result_path+"output_sub_2.h5ad")
