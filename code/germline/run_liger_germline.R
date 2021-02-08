library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
part=args[1]
result_path <- paste0("./results/germline/",part,"/liger/")
dir.create(result_path)
ifile=paste0("./data/germline/",part,"/",part,".rds")
b.combined <- readRDS(ifile)
VariableFeatures(b.combined) <- rownames(b.combined)

b.combined <- ScaleData(b.combined, split.by = "batch", do.center = FALSE)
b.combined <- RunOptimizeALS(b.combined, k = 16, lambda = 2, split.by = "batch")
b.combined <- RunQuantileAlignSNF(b.combined, split.by = "batch")
#b.combined <- RunQuantileNorm(b.combined,split.by = "batch")

b.combined<-RenameAssays(b.combined, RNA = 'integrated')

saveRDS(b.combined, file = paste0(result_path,"integrated_data.rds"))

download.file("")
