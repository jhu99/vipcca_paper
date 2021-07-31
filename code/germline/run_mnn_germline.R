rm(list = ls())
library(Seurat)
#library(SeuratData)
library(SingleCellExperiment)
suppressMessages(library(Matrix))
library(batchelor)
set.seed(42)

args = commandArgs(trailingOnly=TRUE)
part=args[1]
result_path <- paste0("./results/germline/",part,"/mnn/")
dir.create(result_path)
ifile=paste0("./data/germline/",part,"/",part,".rds")
integrated <- readRDS(ifile)
VariableFeatures(integrated) <- rownames(integrated)
genes.use <- rownames(integrated)

Idents(object = integrated) <- "batch"
b1 <- integrated[genes.use, WhichCells(object = integrated, idents = "Guo")]
b2 <- integrated[genes.use, WhichCells(object = integrated, idents = "Li")]

mnnc.results <- mnnCorrect(
	as.SingleCellExperiment(b1),
	as.SingleCellExperiment(b2),
	k=16
)

mdata<-rbind(b1@meta.data,b2@meta.data)
saveRDS(mnnc.results,paste0("./results/germline/",part,"/mnn/mnnc_results.rds"))

mnnc.results<-readRDS(paste0("./results/germline/",part,"/mnn/mnnc_results.rds"))
obj<-CreateSeuratObject(assay(mnnc.results),assay = "integrated",meta.data = mdata)
VariableFeatures(obj)<-rownames(obj)
obj<-ScaleData(object = obj)
obj<-RunPCA(object = obj,npcs = 16,verbose = FALSE)
saveRDS(object = obj, file = paste0("./results/germline/",part,"/mnn/integrated_data.rds"))
