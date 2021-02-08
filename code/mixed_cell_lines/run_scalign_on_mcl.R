library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(SeuratData)
library(PMA)

rm(list = ls())
result_path <- "./results/mixed_cell_lines/real/scalign/"
dir.create(result_path,showWarnings = FALSE)
infile_b1 <- "./data/mixed_cell_lines/293t_filtered_matrices_mex/hg19/"
infile_b2 <- "./data/mixed_cell_lines/jurkat_filtered_matrices_mex/hg19/"
infile_b3 <- "./data/mixed_cell_lines/mixed_filtered_matrices_mex/hg19/"
ann <- read.csv("./data/mixed_cell_lines/293t_jurkat_cluster.txt", header = FALSE)
genelist<- read.csv("./data/mixed_cell_lines/genelist.txt",header = FALSE,stringsAsFactors = FALSE)$V1

b1.data <- Read10X(infile_b1,gene.column = 2)
b2.data <- Read10X(infile_b2,gene.column = 2)
b3.data <- Read10X(infile_b3,gene.column = 2)
b1 <- CreateSeuratObject(counts = b1.data)
b2 <- CreateSeuratObject(counts = b2.data)
b3 <- CreateSeuratObject(counts = b3.data)

b1[["X_batch"]]="293t"
b1[["celltype"]]="293t"
b2[["X_batch"]]="jurkat"
b2[["celltype"]]="jurkat"
b3[["X_batch"]]="mixed"
b3[["celltype"]]=tail(ann, n= 3388)$V1
b1<-RenameCells(b1,add.cell.id = "293t")
b2<-RenameCells(b2,add.cell.id = "jurkat")
b3<-RenameCells(b3,add.cell.id = "mixed")

dataset.list <- list(b1,b2,b3)
for (i in 1:length(dataset.list)) {
	obj <- NormalizeData(dataset.list[[i]], verbose = FALSE, scale.factor = 1e6)
	obj <- ScaleData(obj, do.scale=T, do.center=T, display.progress=T)
	VariableFeatures(obj)<-genelist
	dataset.list[[i]]<-obj
}


## Extract common set of genes across all datasets.
genes.use = Reduce(intersect, lapply(dataset.list, function(seurat.obj) VariableFeatures(seurat.obj)))

## Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
dataset.sce.list = lapply(dataset.list,
                           function(seurat.obj){
                             SingleCellExperiment(assays = list(counts = seurat.obj@assays$RNA@counts[genes.use,],
                                                                logcounts = seurat.obj@assays$RNA@data[genes.use,],
                                                                scale.data = seurat.obj@assays$RNA@scale.data[genes.use,]),
                                                  colData = seurat.obj@meta.data)
                            })
## Create the combined scAlign object from list of SCE(s).
scAlignSce = scAlignCreateObject(sce.objects = dataset.sce.list,
																			genes.use = genes.use,
																			cca.reduce=T,
																			ccs.compute=10,
																			project.name = "scAlign_proj")
scAlignSce = scAlignMulti(scAlignSce,
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
															 run.decoder=FALSE,
															 log.results=TRUE,
															 log.dir=result_path,,
															 device="GPU")
saveRDS(scAlignSce, file = paste0(result_path,"scalign_obj.rds"))

scAlignSce<-readRDS(paste0(result_path,"scalign_obj.rds"))
library(dplyr)
scAlignSeuratObj = as.Seurat(scAlignSce, counts="counts", scale.data="scale.data")
VariableFeatures(scAlignSeuratObj)<-rownames(scAlignSeuratObj)
scAlignSeuratObj <- RunUMAP(scAlignSeuratObj, reduction = "ALIGNED.MultiCCA", dims = 1:16)
saveRDS(data.integrated, file = paste0(result_path,"scalign_seuratobj.rds"))

