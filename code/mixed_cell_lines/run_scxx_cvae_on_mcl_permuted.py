import scxx as sc
import scxx.preprocessing as pp
import scxx.plotting as pl
import sys

mode=str(sys.argv[1])
lambda_regulizer=float(sys.argv[2])
batch_input_size=int(sys.argv[3])
batch_input_size2=int(sys.argv[4])

test_result_path="./results/mixed_cell_lines/permuted/vipcca/%s_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))
import pandas as pd
genelist = pd.read_csv("./data/mixed_cell_lines/genelist_python.txt",header=None,index_col=0).index

r1="./data/mixed_cell_lines/293t.h5ad"
r2="./data/mixed_cell_lines/jurkat.h5ad"
r3="./data/mixed_cell_lines/mixed_permuted.h5ad"

adata_b1 = pp.read_sc_data(r1, batch_name="293t")
adata_b1 = adata_b1[:,genelist]
adata_b2 = pp.read_sc_data(r2, batch_name="jurkat")
adata_b2 = adata_b2[:,genelist]
adata_b3_permuted = pp.read_sc_data(r3, batch_name="permuted")
adata_b3_permuted = adata_b3_permuted[:,genelist]

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=sc.SCXX(
							datasets=[adata_b1,adata_b2,adata_b3_permuted],
							res_path=test_result_path,
							mode=mode,
							split_by="_batch",
							patience_es=50,
							patience_lr=20,
							mt_ratio=0.8,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							batch_input_size2=batch_input_size2,
							# keep_order=True,
							#model_file="model.h5"
							)
adata_transform=handle.fit_transform()
adata_transform.obs.to_csv("./data/mixed_cell_lines/annotation_permuted.txt",sep="\t")
import numpy as np
# 
pl.plotPrediction2(adata_transform.raw.X,adata_transform.X,result_path=test_result_path,rnum=2000,lim=20)
pl.run_embedding(adata_transform,path=test_result_path,method="umap")
pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="_batch",legend_loc="right margin")
pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="celltype",legend_loc="on data")

