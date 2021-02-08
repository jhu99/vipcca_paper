# run seurat on 1M
rm(list = ls())
library(Seurat)
library(liger)
library(SeuratWrappers)
library(SingleCellExperiment)
library(batchelor)
library(Matrix)
result_path="./results/1M_mouse_brain/mnn/"

dir.create(result_path,showWarnings = F)
fns<-dir("./data/scanorama_data/data/mouse_brain",pattern = "data_subset",recursive = T,full.names = T)
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
	# adata<- FindVariableFeatures(adata,selection.method = "vst",nfeatures=2000)
	VariableFeatures(adata)<- rownames(adata)
	datasets<-c(datasets,adata)
}

mnnc.results <- mnnCorrect(
	as.SingleCellExperiment(datasets[[1]]),
	as.SingleCellExperiment(datasets[[2]]),
	as.SingleCellExperiment(datasets[[3]]),
	as.SingleCellExperiment(datasets[[4]]),
	as.SingleCellExperiment(datasets[[5]]),
	as.SingleCellExperiment(datasets[[6]]),
	as.SingleCellExperiment(datasets[[7]]),
	as.SingleCellExperiment(datasets[[8]]),
	as.SingleCellExperiment(datasets[[9]]),
	as.SingleCellExperiment(datasets[[10]]),
	k=16
)
saveRDS(data.combined, file = paste0(result_path,"integrated.rds"))


