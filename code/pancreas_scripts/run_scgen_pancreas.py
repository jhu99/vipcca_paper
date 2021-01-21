import scgen
import scanpy as sc 

test_input_file="./data/panc8/panc8.h5ad"
train = sc.read(test_input_file)
sc.pp.normalize_per_cell(train)
sc.pp.log1p(train)

train.obs["cell_type"] = train.obs["celltype"].tolist()

network = scgen.VAEArith(x_dimension= train.shape[1], model_path="./results/pancreas/scgen/" )
network.train(train_data=train, n_epochs=100)
corrected_adata =  scgen.batch_removal(network, train, batch_key="dataset", cell_label_key="cell_type")
corrected_adata.write("./results/pancreas/scgen/output.h5ad")
