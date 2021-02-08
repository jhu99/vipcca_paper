library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(SeuratData)
library(PMA)

# 
result_path="./results/1M_mouse_brain/scalign/"
# 
dir.create(result_path, showWarning=FALSE)
fns<-dir("./data/scanorama_data/data/mouse_brain",pattern = "data_subset.h5ad",recursive = T,full.names = T)
bnames<-dir("./data/scanorama_data/data/mouse_brain/dropviz/")
bnames<-c(bnames,"spinalCord")
datasets<-list()
i=1
for(fn in fns){
	print(fn)
	adata<-ReadH5AD(fn)
	adata[["batch"]]=bnames[i]
	i=i+1
	adata <- NormalizeData(adata, verbose = FALSE)
	VariableFeatures(adata)<- rownames(adata)
	datasets<-c(datasets,adata)
}


for (i in 1:length(datasets)) {
	obj <- NormalizeData(datasets[[i]], verbose = FALSE, scale.factor = 1e6)
	obj <- ScaleData(obj, do.scale=T, do.center=T, display.progress=T)
	#datasets[[i]] <- FindVariableFeatures(datasets[[i]], selection.method = "vst",
	#																					 nfeatures = 2000, verbose = FALSE)
	VariableFeatures(obj)<-rownames(obj)
	datasets[[i]]<-obj
}

## Extract common set of genes across all datasets.
genes.use = Reduce(intersect, lapply(datasets, function(seurat.obj) VariableFeatures(seurat.obj)))

## Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
pancreas.sce.list = lapply(datasets,
                           function(seurat.obj){
                             SingleCellExperiment(assays = list(counts = seurat.obj@assays$RNA@counts[genes.use,],
                                                                logcounts = seurat.obj@assays$RNA@data[genes.use,],
                                                                scale.data = seurat.obj@assays$RNA@scale.data[genes.use,]),
                                                  colData = seurat.obj@meta.data)
                            })
## Create the combined scAlign object from list of SCE(s).
scAlignPancreas = scAlignCreateObject(sce.objects = pancreas.sce.list,
																			genes.use = genes.use,
																			cca.reduce=T,
																			ccs.compute=10,
																			project.name = "scAlign_Pancreatic_Islet")
scAlignPancreas = scAlignMulti(scAlignPancreas,
															 options=scAlignOptions(steps=1500,
															 											 log.every=500,
															 											 batch.size=300,
															 											 perplexity=30,
															 											 norm=TRUE,
															 											 batch.norm.layer=FALSE,
															 											 architecture="large",  ## 3 layer neural network
															 											 num.dim=16),            ## Number of latent dimensions
															 encoder.data="MultiCCA",
															 decoder.data="scale.data",
															 supervised='none',
															 run.encoder=TRUE,
															 run.decoder=TRUE,
															 log.results=TRUE,
															 log.dir=result_path,,
															 device="GPU")
saveRDS(scAlignPancreas, file = paste0(result_path,"scalign_pancreas.rds"))

scAlignPancreas=readRDS(paste0(result_path,"scalign_pancreas.rds"))
library(dplyr)
scAlignSeuratObj = as.Seurat(scAlignPancreas, counts="counts", scale.data="scale.data")
pancreas.integrated <- RunUMAP(scAlignSeuratObj, reduction = "ALIGNED.MultiCCA", dims = 1:16)
saveRDS(pancreas.integrated, file = paste0(result_path,"scalign_pancreas_integrated.rds"))

