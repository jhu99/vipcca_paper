# run Seurat V3 on bipolar
library(Seurat)
library(dplyr)
library(kBET)
library(ggplot2)

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
part=args[1]
result_path <- paste0("./results/germline/",part,"/uncorrected/")
dir.create(result_path)
ifile=paste0("./data/germline/",part,"/",part,".rds")
integrated <- readRDS(ifile)
VariableFeatures(integrated) <- rownames(integrated)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 16, verbose = FALSE)

integrated<-RenameAssays(integrated, RNA = 'integrated')
saveRDS(integrated,file =  paste0(result_path,"integrated_data.rds"))

