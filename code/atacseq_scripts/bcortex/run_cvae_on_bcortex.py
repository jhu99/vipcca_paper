import scxx as sc
import scxx.preprocessing as pp
import scxx.plotting as pl
import sys
import numpy as np
import pandas as pd
from keras import optimizers

mode=str(sys.argv[1])
lambda_regulizer=float(sys.argv[2])
batch_input_size=int(sys.argv[3])
batch_input_size2=int(sys.argv[4])
test_input_file1="./data/atacseq/snareseq/12A/ExpressionMatrix_cDNA/"
test_input_file2="./data/atacseq/snareseq/12A/GeneActivity_chromatin/"
test_result_path="./results/atacseq/12A/scxx_%s/scmb_%s_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))

ann = pd.read_csv("./data/atacseq/snareseq/12A/annotation.csv",sep=" ")
adata1 = pp.read_sc_data(test_input_file1,batch_name="ATAC",fmt="10x_mtx")
adata1.obs=adata1.obs.join(ann)
adata2 = pp.read_sc_data(test_input_file2,batch_name="RNA",fmt="10x_mtx")
adata2.obs=adata2.obs.join(ann)

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=sc.SCXX(
							datasets=[adata1, adata2],
							res_path=test_result_path,
							mode=mode,
							split_by="_batch",
							patience_es=100,
							patience_lr=50,
							mt_ratio=0.8,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							batch_input_size2=batch_input_size2,
							hidden_layers=[128,64,32,16],
							index_unique="-",
							#optimizer=optimizers.RMSprop(lr=0.01),
							model_file="model.h5"
							)
adata_transform=handle.fit_transform()
adata_transform.obs.to_csv("./data/atacseq/snareseq/12A/annotation_vipcca.csv")
#
import numpy as np
np.savetxt("./data/atacseq/snareseq/12A/genelist_python.txt",adata_transform.var_names.values,fmt="%s")

# pl.plotPrediction2(adata_transform.raw.X,adata_transform.X,result_path=test_result_path,lim=20)
# pl.run_embedding(adata_transform,path=test_result_path,method="umap")
# pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="_batch",legend_loc="right margin")
#pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="celltype",legend_loc="on data")
# pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="louvain",legend_loc="right margin")
