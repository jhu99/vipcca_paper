import scanpy as sc
import desc
import os
import sys

part=str(sys.argv[1])
test_result_path="./results/germline/"+part+"/desc/"
rfile="./data/germline/"+part+"/"+part+".h5ad"
if not os.path.exists(test_result_path):
	os.makedirs(test_result_path)

adata=sc.read_h5ad(rfile)

sc.pp.log1p(adata)
adataw=adata.copy()
sc.pp.scale(adata,max_value=6)

adata=desc.train(adata,
        dims=[adata.shape[1],32,16],
        tol=0.005,
        n_neighbors=10,
        batch_size=50,
        louvain_resolution=0.45,# not necessarily a list, you can only set one value, like, louvain_resolution=1.0
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
