import scxx as sc
import scxx.preprocessing as pp
import scxx.plotting as pl
import scxx.network as network
import sys
import numpy as np
import pandas as pd
from keras import optimizers

mode=str(sys.argv[1])
lambda_regulizer=float(sys.argv[2])
batch_input_size=int(sys.argv[3])
batch_input_size2=int(sys.argv[4])
test_input_file1="./data/atacseq/10x_pbmc/ExpressMat/"
test_input_file2="./data/atacseq/10x_pbmc/geneActMat/"
test_result_path="./results/atacseq/10x_pbmc/scxx_%s/scmb_%s_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))

# ifile="./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_3_128_10/output.h5ad"
# genelist = pd.read_csv("./data/atacseq/10x_pbmc/genelist.txt", header=None, index_col=0)
# adata1 = pp.read_sc_data(test_input_file1,batch_name="RNA",fmt="10x_mtx")
# adata1 = adata1[:,genelist.index]
# adata2 = pp.read_sc_data(test_input_file2,batch_name="ATAC",fmt="10x_mtx")
# adata2 = adata2[:,genelist.index]
ann = pd.read_csv("./data/atacseq/10x_pbmc/annotation.csv",index_col=0)
adata1 = pp.read_sc_data("./data/atacseq/10x_pbmc/ExpressMat/rna.h5ad")
adata2 = pp.read_sc_data("./data/atacseq/10x_pbmc/geneActMat/geneact.h5ad")
# adata1.obs = adata1.obs.join(ann)
# adata2.obs = adata2.obs.join(ann)
# adata1.write("./data/atacseq/10x_pbmc/ExpressMat/rna.h5ad")
# adata2.write("./data/atacseq/10x_pbmc/geneActMat/geneact.h5ad")

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=sc.SCXX(
							datasets=[adata1, adata2],
							res_path=test_result_path,
							mode=mode,
							split_by="_batch",
							patience_es=20,
							patience_lr=10,
							mt_ratio=0.8,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							batch_input_size2=batch_input_size2,
							hidden_layers=[128,64,32,16],
							hvg=False,
							model_file="model.h5",
							save=False
							)
adata_transform=handle.fit_transform()
ann_atac=network.findNeighbors(adata_transform)
adata_transform.write(test_result_path+"output.h5ad")
ann_atac.to_csv(test_result_path+"ann_atac.csv")
adata_transform.obs.to_csv(test_result_path+"annotation.csv")



