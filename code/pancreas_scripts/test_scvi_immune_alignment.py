
n_epochs_all = None
save_path = 'data/immune_alignment/'
show_plot = True

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

import numpy as np
import numpy.random as random
import pandas as pd
import scanpy as sc

use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE


from umap import UMAP
from scvi.dataset.csv import CsvDataset
stimulated = CsvDataset(filename='immune_stimulated_expression_matrix.txt',
                        save_path=save_path,sep='\t', new_n_genes=35635)
control = CsvDataset(filename='immune_control_expression_matrix.txt',
                     save_path=save_path, sep='\t', new_n_genes=35635)
                     
all_dataset = GeneExpressionDataset()
all_dataset.populate_from_datasets([control, stimulated])
vae = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels,
          n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')

trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0,use_cuda=use_cuda)
n_epochs = 100 if n_epochs_all is None else n_epochs_all
trainer.train(n_epochs=n_epochs)


