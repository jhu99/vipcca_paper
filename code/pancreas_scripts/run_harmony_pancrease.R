rm(list = ls())
library(SeuratData)
library(Seurat)
library(harmony)
set.seed(0)
data("panc8")
result_path="./results/pancreas/harmony/"
dir.create(result_path, showWarning=FALSE)
pancreas.list <- SplitObject(panc8, split.by = "dataset")
genelist<- read.csv("./data/panc8/genelist.txt",header = FALSE,stringsAsFactors = FALSE)$V1
for (i in 1:length(pancreas.list)) {
	obj <- NormalizeData(pancreas.list[[i]], verbose = FALSE, scale.factor = 1e6)
	VariableFeatures(object = obj)<-genelist
	pancreas.list[[i]] <- obj
}
adata<-merge(pancreas.list[[1]],y=pancreas.list[2:length(pancreas.list)], merge.data=TRUE)
VariableFeatures(adata)<-genelist
remove(panc8,pancreas.list)
adata<-ScaleData(adata)
adata<-RunPCA(adata,npcs = 16)

adata.harmony<-RunHarmony(adata,"tech")
saveRDS(adata.harmony, paste0(result_path,"integrated_data.rds"))
library(ggplot2)
library(cowplot)
adata.harmony<-RunUMAP(adata.harmony,reduction = "harmony",dims = 1:16)
DimPlot(adata.harmony,reduction = "umap",group.by = "celltype")
DimPlot(adata.harmony,reduction = "umap",group.by = "dataset")
