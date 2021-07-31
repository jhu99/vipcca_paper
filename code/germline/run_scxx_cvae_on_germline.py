import VIPCCA as vp
from VIPCCA import preprocessing as pp
from VIPCCA import plotting as pl
import scanpy
import sys

mode=str(sys.argv[1])
lambda_regulizer=float(sys.argv[2])
batch_input_size=int(sys.argv[3])
batch_input_size2=int(sys.argv[4])
gender=str(sys.argv[5])
repeat=int(sys.argv[6])

# Read and preproces datasets.
test_result_path="./results/germline/%s/VIPCCA/%s_%s_%s_%s_%s/" % (str(sys.argv[5]).lower(),str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),str(sys.argv[6]))
rfile="./data/germline/"+gender+"/"+gender+".h5ad"
adata=pp.read_sc_data(rfile)
scanpy.pp.log1p(adata)

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=vp.VIPCCA(
							adata,
							res_path=test_result_path,
							mode=mode,
							split_by="batch",
							epochs=100,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							batch_input_size2=batch_input_size2,
							# model_file="model.h5"
							) 
adata_transform=handle.fit_transform()
pl.plotPrediction2(adata_transform.raw.X,adata_transform.X,result_path=test_result_path, rnum=80, lim=20)

