# run Seurat V3 on bipolar
library(Seurat)
library(dplyr)
library(kBET)
library(ggplot2)
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
part=args[1]
result_path <- paste0("./results/germline/",part,"/seurat/")
dir.create(result_path)
ifile=paste0("./data/germline/",part,"/",part,".rds")
adata <- readRDS(ifile)
VariableFeatures(adata) <- rownames(adata)

datasets <-	SplitObject(adata,split.by = "batch")
b.anchors <- FindIntegrationAnchors(object.list = datasets, dims = 1:16,anchor.features=rownames(adata),k.filter = 50)
b.integrated <- IntegrateData(anchorset = b.anchors, dims = 1:16)

library(ggplot2)
library(cowplot)
DefaultAssay(b.integrated) <- "integrated"
b.integrated <- ScaleData(b.integrated, verbose = FALSE)
b.integrated <- RunPCA(b.integrated, npcs = 16, verbose = FALSE)

saveRDS(b.integrated,file = paste0(result_path,"integrated_data.rds"))
