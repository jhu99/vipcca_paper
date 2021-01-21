import scbean.model.vipcca as vip
import scbean.tools.utils as tl
import scbean.tools.plotting as pl

# import sys
# mode=str(sys.argv[1])
# lambda_regulizer=float(sys.argv[2])
# batch_input_size=int(sys.argv[3])
# batch_input_size2=int(sys.argv[4])

test_input_file="../data/panc8/panc8.h5ad"
test_result_path="../results/pancreas/scxx_%s_1e4/%s_var2_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))

adata = tl.read_sc_data(test_input_file)
datasets=tl.split_object(adata,by="dataset")
adata_all = tl.preprocessing(datasets)
# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=vip.VIPCCA(
							    adata_all=adata_all,
                  res_path=test_result_path,
                  mode='CVAE',
                  split_by="dataset",
                  patience_es=50,
                  patience_lr=20,
                  lambda_regulizer=5,
                  batch_input_size=128,
                  batch_input_size2=16,
                  epochs=20,
							)
adata_integrate=handle.fit_integrate()
import numpy as np
np.savetxt("../data/panc8/genelist.txt",adata_integrate.var_names.values,fmt="%s")

pl.plotCorrelation(adata_integrate.raw.X,adata_integrate.X,result_path=test_result_path,lim=20)
pl.run_embedding(adata_integrate,path=test_result_path,method="umap")
pl.plotEmbedding(adata_integrate, path=test_result_path, method='umap', group_by="dataset",legend_loc="right margin")
pl.plotEmbedding(adata_integrate, path=test_result_path, method='umap', group_by="celltype",legend_loc="on data")
