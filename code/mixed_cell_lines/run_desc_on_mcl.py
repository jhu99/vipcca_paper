import scanpy as sc
import scxx.preprocessing as pp
import desc
import os

test_result_path="./results/mixed_cell_lines/real/desc/"
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
        num_Cores=1, #for reproducible, only use 1 cpu
        num_Cores_tsne=4,
        save_encoder_weights=False,
        save_encoder_step=3,# save_encoder_weights is False, this parameter is not used
        use_ae_weights=False,
        do_umap=False)
adata.raw=adataw
adata.write(test_result_path+"output.h5ad")
