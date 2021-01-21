library(SeuratData)
library(Seurat)

data("panc8")
result_path="./results/pancreas/seurat/"
dir.create(result_path, showWarning=FALSE)
pancreas.list <- SplitObject(panc8, split.by = "dataset")
genelist<- read.csv("./data/panc8/genelist.txt",header = FALSE,stringsAsFactors = FALSE)$V1
for (i in 1:length(pancreas.list)) {
	obj <- NormalizeData(pancreas.list[[i]], verbose = FALSE, scale.factor = 1e6)
	VariableFeatures(object = obj)<-genelist
	pancreas.list[[i]] <- obj
	#pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst",
	# 																					 nfeatures = 2000, verbose = FALSE)
}

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:16)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:16)

library(ggplot2)
library(cowplot)
# # switch to integrated assay. The variable features of this assay are automatically
# # set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 16, verbose = FALSE)

saveRDS(pancreas.integrated,file = "./results/pancreas/seurat/seurat_v3_integrated_data.rds")

