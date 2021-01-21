import scxx as sc
import scxx.preprocessing as pp
import scxx.plotting as pl
import sys
from keras import optimizers

mode=str(sys.argv[1])
lambda_regulizer=float(sys.argv[2])
batch_input_size=int(sys.argv[3])
test_input_file="./data/panc8/panc8.h5ad"
test_result_path="./results/pancreas/scxx_%s/scmb_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]))

adata = pp.read_sc_data(test_input_file)
datasets=pp.split_object(adata,by="dataset")

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=sc.SCXX(
							rawdata=adata,
							datasets=datasets,
							res_path=test_result_path,
							mode=mode,
							split_by="dataset",
							patience_es=1000,
							patience_lr=50,
							mt_ratio=0.8,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							keep_order=True,
							hidden_layers=[128,64,32,16],
							optimizer=optimizers.RMSprop(lr=0.01),
							#model_file="model.h5"
							)
adata_transform=handle.fit_transform()
# 
pl.plotPrediction2(adata_transform.raw.X,adata_transform.X,result_path=test_result_path,lim=500)
# pl.run_embedding(adata_transform,path=test_result_path,method="umap")
# pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="dataset",legend_loc="right margin")
# pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="celltype",legend_loc="on data")
# pl.plotEmbedding(adata_transform, path=test_result_path, method='umap', group_by="louvain",legend_loc="right margin")
