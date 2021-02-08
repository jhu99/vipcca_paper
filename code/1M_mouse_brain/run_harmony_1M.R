# run seurat on 1M
rm(list = ls())
library(Seurat)
library(harmony)

result_path="./results/1M_mouse_brain/harmony/"
genelist<- read.csv("./data/scanorama_data/data/mouse_brain/genelist_vipcca.txt",header = FALSE,stringsAsFactors = FALSE)$V1
dir.create("./results/1M_mouse_brain/harmony/",showWarnings = F)
fns<-dir("./data/scanorama_data/data/mouse_brain",pattern = "data.h5ad",recursive = T,full.names = T)
bnames<-dir("./data/scanorama_data/data/mouse_brain/dropviz/")
bnames<-c("spinalCord",bnames)
fns
datasets<-list()
i=1
for(fn in fns){
	print(fn)
	adata<-ReadH5AD(fn)
	
	adata[[".batch"]]=i
	i=i+1
	adata <- NormalizeData(adata, verbose = FALSE)
	adata<-adata[genelist,]
	VariableFeatures(adata)<- rownames(adata)
	datasets<-c(datasets,adata)
}

adata<-merge(datasets[[1]],y=datasets[2:length(datasets)], merge.data=TRUE)
adata$.batch<- as.factor(adata$.batch)
remove(datasets)

VariableFeatures(adata)<-rownames(adata)
adata<-ScaleData(adata)
adata<-RunPCA(adata,npcs = 16)

adata.harmony<-RunHarmony(adata,".batch")
adata.harmony <- DietSeurat(adata.harmony, counts = TRUE, data = TRUE, scale.data = FALSE,dimreducs = "harmony")
# saveRDS(adata.harmony,"./results/1M_mouse_brain/harmony/integrated_data.rds")

adata.harmony<-readRDS("./results/1M_mouse_brain/harmony/integrated_data.rds")
adata.harmony$.batch<-as.integer.Array(adata.harmony$.batch)
ann<-adata.harmony@meta.data
ann[ann$.batch==10,1]=0
adata.harmony<-AddMetaData(adata.harmony,metadata = ann)
saveRDS(adata.harmony, paste0(result_path,"integrated_data.rds"))



