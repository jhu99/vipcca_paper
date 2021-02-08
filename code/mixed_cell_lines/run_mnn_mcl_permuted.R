library(Seurat)
library(SeuratData)
library(SingleCellExperiment)
suppressMessages(library(Matrix))
library(batchelor)

rm(list = ls())
set.seed(42)
result_path="./results/mixed_cell_lines/permuted/mnn/"
dir.create(result_path, showWarning=FALSE)
infile_b1 <- "./data/mixed_cell_lines/293t_filtered_matrices_mex/hg19/"
infile_b2 <- "./data/mixed_cell_lines/jurkat_filtered_matrices_mex/hg19/"
infile_b3 <- "./data/mixed_cell_lines/mixed_filtered_matrices_mex_permuted/"
ann <- read.csv("./data/mixed_cell_lines/293t_jurkat_cluster.txt", header = FALSE)
genelist<- read.csv("./data/mixed_cell_lines/genelist_permuted.txt",header = FALSE,stringsAsFactors = FALSE)$V1

b1.data <- Read10X(infile_b1,gene.column = 2)
b2.data <- Read10X(infile_b2,gene.column = 2)
b3.data <- Read10X(infile_b3,gene.column = 1)
b1 <- CreateSeuratObject(counts = b1.data)
b2 <- CreateSeuratObject(counts = b2.data)
b3 <- CreateSeuratObject(counts = b3.data)

b1[["X_batch"]]="293t"
b1[["celltype"]]="293t"
b2[["X_batch"]]="jurkat"
b2[["celltype"]]="jurkat"
b3[["X_batch"]]="mixed_permuted"
b3[["celltype"]]=tail(ann, n= 3388)$V1
b1<-RenameCells(b1,add.cell.id = "293t")
b2<-RenameCells(b2,add.cell.id = "jurkat")
b3<-RenameCells(b3,add.cell.id = "mixed_permuted")

b1 <- NormalizeData(b1,verbose = TRUE,scale.factor = 1e6)
b2 <- NormalizeData(b2,verbose = TRUE,scale.factor = 1e6)
b3 <- NormalizeData(b3,verbose = TRUE,scale.factor = 1e6)
b1 <- b1[genelist,]
b2 <- b2[genelist,]
b3 <- b3[genelist,]

mnnc.results <- mnnCorrect(
	as.SingleCellExperiment(b1),
	as.SingleCellExperiment(b2),
	as.SingleCellExperiment(b3),
	k=16
)
mdata<-rbind(b1@meta.data,b2@meta.data,b3@meta.data)

saveRDS(mnnc.results,"../results/mixed_cell_lines/real/mnn/mnnc_results.rds")

mnnc.results<-readRDS("./results/mixed_cell_lines/real/mnn/mnnc_results.rds")
obj<-CreateSeuratObject(assay(mnnc.results),assay = "integrated",meta.data = mdata)
VariableFeatures(obj)<-rownames(obj)
obj<-ScaleData(object = obj)
obj<-RunPCA(object = obj,npcs = 16,verbose = FALSE)
saveRDS(object = obj, file = "./results/mixed_cell_lines/real/mnn/mnn_integrated_data.rds")
