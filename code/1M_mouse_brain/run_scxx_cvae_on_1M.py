import scxx as sc
import scxx.preprocessing as pp
import scxx.plotting as pl
import sys
import numpy as np
#from process import laod_names

mode=str(sys.argv[1])
lambda_regulizer=float(sys.argv[2])
batch_input_size=int(sys.argv[3])
batch_input_size2=int(sys.argv[4])

test_result_path="./results/1M_mouse_brain/VIPCCA/%s_%s_%s_%s/" % (str(sys.argv[1]).lower(),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))

data_names = [
    'data/mouse_brain/nuclei',
    'data/mouse_brain/dropviz/Cerebellum_ALT',
    'data/mouse_brain/dropviz/Cortex_noRep5_FRONTALonly',
    'data/mouse_brain/dropviz/Cortex_noRep5_POSTERIORonly',
    'data/mouse_brain/dropviz/EntoPeduncular',
    'data/mouse_brain/dropviz/GlobusPallidus',
    'data/mouse_brain/dropviz/Hippocampus',
    'data/mouse_brain/dropviz/Striatum',
    'data/mouse_brain/dropviz/SubstantiaNigra',
    'data/mouse_brain/dropviz/Thalamus',
]

datasets=[]
for i in range(len(data_names)):
	name=data_names[i]
	ann = pp.read_sc_data("./data/scanorama_data/"+name+"/data.h5ad",batch_name=str(i))
	datasets.append(ann)

# Run scxx. Here, rawdata is optional. It will be generated from datasets if not provided.
handle=sc.SCXX(
							datasets=datasets,
							res_path=test_result_path,
							mode=mode,
							split_by="_batch",
							patience_es=10,
							patience_lr=5,
							mt_ratio=0.8,
							lambda_regulizer=lambda_regulizer,
							batch_input_size=batch_input_size,
							batch_input_size2=batch_input_size2,
							index_unique=None,
							# keep_order=True,
							# model_file="model.h5"
							)
adata_transform=handle.fit_transform()
import numpy as np
np.savetxt("./data/scanorama_data/data/mouse_brain/genelist_vipcca.txt",adata_transform.var_names.values,fmt="%s")

# 
# sampling and visualization
# pl.plotPrediction2(adata_transform.raw.X,adata_transform.X,result_path=test_result_path,rnum=2000,lim=20)
cell_ind=np.random.choice(adata_transform.shape[0],20000,replace=False)
adata_transform_sub=adata_transform[adata_transform.obs.index[cell_ind],:]
adata_transform_sub.write("./results/1M_mouse_brain/VIPCCA/cvae_2_64_16/output_sub.h5ad")
adata_transform_sub.obs.to_csv("./results/1M_mouse_brain/VIPCCA/cvae_2_64_16/annotation.csv")

