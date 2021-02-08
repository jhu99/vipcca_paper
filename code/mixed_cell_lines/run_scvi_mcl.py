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

test_result_path="./results/mixed_cell_lines/real/scvi/"
import pandas as pd
genelist = pd.read_csv("./data/mixed_cell_lines/genelist_python.txt",header=None,index_col=0).index

if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)

r1="./data/mixed_cell_lines/293t.h5ad"
r2="./data/mixed_cell_lines/jurkat.h5ad"
r4="./data/mixed_cell_lines/mixed.h5ad"

adata_b1 = pp.read_sc_data(r1, batch_name="293t")
adata_b2 = pp.read_sc_data(r2, batch_name="jurkat")
adata_b4 = pp.read_sc_data(r4, batch_name="mixed")
adata=adata_b1.concatenate(adata_b2,adata_b4,index_unique="-")
adata=adata[:,genelist]

#adata.obs["_batch"]=adata.obs._batch.astype("category").cat.codes
scviDataset = AnnDatasetFromAnnData(adata, batch_label="batch")
model_file=test_result_path + "model.pkl"


use_cuda=True
use_batches=True
n_epochs = 100
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
adata_s=AnnData(X=cell_freq,obs=adata.obs,var=adata.var)
cell_freq[cell_freq<1e-6]=0
adata_s.raw=AnnData(X=csr_matrix(cell_freq),obs=adata.obs,var=adata.var)
adata_s.obsm['X_scvi']=latent
adata_s.obs.to_csv("./data/mixed_cell_lines/annotation_scvi.txt",sep="\t")
adata_s.write(test_result_path+"output.h5ad")
