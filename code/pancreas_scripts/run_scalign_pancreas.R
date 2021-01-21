library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(SeuratData)
library(PMA)

data("panc8")
# 
result_path="./results/pancreas/scalign/"
# 
dir.create(result_path, showWarning=FALSE)
genelist<- read.csv("./data/panc8/genelist.txt",header = FALSE,stringsAsFactors = FALSE)$V1
pancreas.list <- SplitObject(panc8, split.by = "dataset")
for (i in 1:length(pancreas.list)) {
	obj <- NormalizeData(pancreas.list[[i]], verbose = FALSE, scale.factor = 1e6)
	obj <- ScaleData(obj, do.scale=T, do.center=T, display.progress=T)
	#pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst",
	#																					 nfeatures = 2000, verbose = FALSE)
	VariableFeatures(obj)<-genelist
	pancreas.list[[i]]<-obj
}

## Extract common set of genes across all datasets.
genes.use = Reduce(intersect, lapply(pancreas.list, function(seurat.obj) VariableFeatures(seurat.obj)))

## Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
pancreas.sce.list = lapply(pancreas.list,
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
															 options=scAlignOptions(steps=15000,
															 											 log.every=5000,
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

