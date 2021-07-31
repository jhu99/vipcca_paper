# run Seurat V3 on bipolar
library(Seurat)
library(dplyr)
library(kBET)
library(ggplot2)
library(harmony)

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
part=args[1]
result_path <- paste0("./results/germline/",part,"/harmony/")
dir.create(result_path)
ifile=paste0("./data/germline/",part,"/",part,".rds")
adata <- readRDS(ifile)
VariableFeatures(adata) <- rownames(adata)

adata<-ScaleData(adata)
adata<-RunPCA(adata,npcs = 16)
adata.harmony<-RunHarmony(adata,"batch")
adata.harmony<-RenameAssays(adata.harmony, RNA = 'integrated')

saveRDS(adata.harmony, paste0(result_path,"integrated_data.rds"))

