rm(list = ls())
.libPaths()
library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(SeuratData)
library(PMA)

args = commandArgs(trailingOnly=TRUE)
part=args[1]
result_path <- paste0("./results/germline/",part,"/scalign/")
dir.create(result_path)
ifile=paste0("./data/germline/",part,"/",part,".rds")
adata <- readRDS(ifile)
VariableFeatures(adata) <- rownames(adata)
reference.list <- SplitObject(adata, split.by = "batch")

for (i in 1:length(reference.list)) {
	# reference.list[[i]] <- NormalizeData(reference.list[[i]], verbose = FALSE)
	reference.list[[i]] <- ScaleData(reference.list[[i]], do.scale=T, do.center=T, display.progress=T)
	VariableFeatures(reference.list[[i]]) <- rownames(reference.list[[i]])
}

## Extract common set of genes across all datasets.
genes.use = Reduce(intersect, lapply(reference.list, function(seurat.obj) VariableFeatures(seurat.obj)))

## Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
reference.sce.list = lapply(reference.list,
                           function(seurat.obj){
                             SingleCellExperiment(assays = list(counts = seurat.obj@assays$RNA@counts[genes.use,],
                                                                logcounts = seurat.obj@assays$RNA@data[genes.use,],
                                                                scale.data = seurat.obj@assays$RNA@scale.data[genes.use,]),
                                                  colData = seurat.obj@meta.data)
                            })
## Create the combined scAlign object from list of SCE(s).
scAlignSce = scAlignCreateObject(sce.objects = reference.sce.list,
																			genes.use = genes.use,
																			cca.reduce=T,
																			ccs.compute=10,
																			project.name = "scAlign_bipolar_Islet")
scAlignSce = scAlignMulti(scAlignSce,
															 options=scAlignOptions(steps=15000,
															 											 log.every=5000,
															 											 batch.size=50,
															 											 perplexity=30,
															 											 norm=TRUE,
															 											 batch.norm.layer=FALSE,
															 											 architecture="large",  ## 3 layer neural network
															 											 num.dim=16),            ## Number of latent dimensions
															 encoder.data="CCA",
															 decoder.data="scale.data",
															 supervised='none',
															 run.encoder=TRUE,
															 run.decoder=TRUE,
															 log.results=TRUE,
															 log.dir=result_path,
															 device="GPU")
saveRDS(scAlignSce, file = paste0(result_path,"scalign_integrated.rds"))

scAlignSce<-readRDS(paste0(result_path,"scalign_integrated.rds"))
library(dplyr)
scAlignSeuratObj = as.Seurat(scAlignSce, counts="counts", scale.data="scale.data")
scAlignSeuratObj[["integrated"]] <- scAlignSeuratObj[["RNA"]]
DefaultAssay(scAlignSeuratObj)<-"integrated"
scAlignSeuratObj<-ScaleData(scAlignSeuratObj)

saveRDS(scAlignSeuratObj, file = paste0(result_path,"integrated_data.rds"))

