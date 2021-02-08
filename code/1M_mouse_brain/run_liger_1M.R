# run seurat on 1M
rm(list = ls())
library(Seurat)
library(liger)
library(SeuratWrappers)
result_path="./results/1M_mouse_brain/liger/"

dir.create(result_path,showWarnings = F)
fns<-dir("./data/scanorama_data/data/mouse_brain",pattern = "data_uniq",recursive = T,full.names = T)
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
	adata<- FindVariableFeatures(adata,selection.method = "vst",nfeatures=2000)
	#VariableFeatures(adata)<- rownames(adata)
	datasets<-c(datasets,adata)
}

data.combined <- merge(datasets[[1]],y=datasets[2:10],add.cell.ids = 1:10)
data.combined <- NormalizeData(data.combined,scale.factor = 1e6)
VariableFeatures(data.combined) <- rownames(data.combined)
data.combined <- ScaleData(data.combined, split.by = "batch", do.center = FALSE)
data.combined <- RunOptimizeALS(data.combined, k = 16, lambda = 5, split.by = "batch")
data.combined <- RunQuantileAlignSNF(data.combined, split.by = "batch", resolution = 0.3)
data.combined <- RunUMAP(data.combined, dims = 1:ncol(data.combined[["iNMF"]]), reduction = "iNMF")

saveRDS(data.combined, file = paste0(result_path,"integrated.rds"))


