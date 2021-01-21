library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)

data("panc8")
result_path="./results/pancreas/liger/"
dir.create(result_path, showWarning=FALSE)
panc8 <- NormalizeData(panc8, scale.factor=1e6)
#panc8 <- FindVariableFeatures(panc8)
genelist<- read.csv("./data/panc8/genelist.txt",header = FALSE,stringsAsFactors = FALSE)$V1
VariableFeatures(panc8)<-genelist
panc8 <- ScaleData(panc8, split.by = "dataset", do.center = FALSE)
panc8 <- RunOptimizeALS(panc8, k = 16, lambda = 5, split.by = "dataset")
panc8 <- RunQuantileAlignSNF(panc8, split.by = "dataset", resolution = 0.5)
panc8 <- RunUMAP(panc8, dims = 1:ncol(panc8[["iNMF"]]), reduction = "iNMF")
saveRDS(panc8, file = paste0(result_path,"liger_integrated_data.rds"))
