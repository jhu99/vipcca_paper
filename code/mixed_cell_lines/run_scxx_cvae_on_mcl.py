import scxx as sc
import scxx.preprocessing as pp
import scxx.plotting as pl
import sys

# "CVAE"
mode=str(sys.argv[1])
# 5
lambda_regulizer=int(sys.argv[2])
# 128
batch_input_size=int(sys.argv[3])
# 16
batch_input_size2=int(sys.argv[4])

test_result_path="./results/mixed_cell_lines/real/vipcca/%s_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))

r1="./data/mixed_cell_lines/293t.h5ad"
r2="./data/mixed_cell_lines/jurkat.h5ad"
r4="./data/mixed_cell_lines/mixed.h5ad"

adata_b1 = pp.read_sc_data(r1, batch_name="293t")
adata_b2 = pp.read_sc_data(r2, batch_name="jurkat")
adata_b4 = pp.read_sc_data(r4, batch_name="mixed")

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=sc.SCXX(
							datasets=[adata_b1,adata_b2,adata_b4],
							res_path=test_result_path,
							mode=mode,
							split_by="_batch",
							patience_es=50,
							patience_lr=20,
							mt_ratio=0.8,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							batch_input_size2=batch_input_size2,
							index_unique="-",
							# keep_order=True,
							model_file="model.h5"
							)
adata_transform=handle.fit_transform()
adata_transform.obs.to_csv("./data/mixed_cell_lines/annotation.txt",sep="\t")
import numpy as np
np.savetxt("./data/mixed_cell_lines/genelist_python.txt",adata_transform.var_names.values,fmt="%s")

# 
pl.plotPrediction2(adata_transform.raw.X,adata_transform.X,result_path=test_result_path,rnum=2000,lim=20)
pl.run_embedding(adata_transform,path=test_result_path,method="umap")
pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="_batch",legend_loc="right margin")
pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="celltype",legend_loc="on data")

