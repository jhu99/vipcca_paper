# run seurat on 1M
rm(list = ls())
library(Seurat)
dir.create("./results/1M_mouse_brain/seurat/",showWarnings = F)
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
integration.anchors<-FindIntegrationAnchors(object.list = datasets,dims = 16)
integrated <- IntegrateData(anchorset = integration.anchors, dims = 1:16)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 16, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:16)
saveRDS(integrated,"./results/1M_mouse_brain/seurat/integrated.rds")

