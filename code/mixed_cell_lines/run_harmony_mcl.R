# run Seurat V3 on bipolar
library(Seurat)
library(dplyr)
library(kBET)
library(ggplot2)
library(harmony)

rm(list = ls())
result_path <- "./results/mixed_cell_lines/real/harmony/"
dir.create(result_path,showWarnings = FALSE)
infile_b1 <- "./data/mixed_cell_lines/293t_filtered_matrices_mex/hg19/"
infile_b2 <- "./data/mixed_cell_lines/jurkat_filtered_matrices_mex/hg19/"
infile_b3 <- "./data/mixed_cell_lines/mixed_filtered_matrices_mex/hg19/"
ann <- read.csv("./data/mixed_cell_lines/293t_jurkat_cluster.txt", header = FALSE)
genelist<- read.csv("./data/mixed_cell_lines/genelist.txt",header = FALSE,stringsAsFactors = FALSE)$V1

b1.data <- Read10X(infile_b1,gene.column = 2)
b2.data <- Read10X(infile_b2,gene.column = 2)
b3.data <- Read10X(infile_b3,gene.column = 2)
b1 <- CreateSeuratObject(counts = b1.data)
b2 <- CreateSeuratObject(counts = b2.data)
b3 <- CreateSeuratObject(counts = b3.data)

b1[["X_batch"]]="293t"
b1[["celltype"]]="293t"
b2[["X_batch"]]="jurkat"
b2[["celltype"]]="jurkat"
b3[["X_batch"]]="mixed"
b3[["celltype"]]=tail(ann, n= 3388)$V1
b1<-RenameCells(b1,add.cell.id = "293t")
b2<-RenameCells(b2,add.cell.id = "jurkat")
b3<-RenameCells(b3,add.cell.id = "mixed")

b1 <- NormalizeData(b1,verbose = TRUE,scale.factor = 1e6)
b2 <- NormalizeData(b2,verbose = TRUE,scale.factor = 1e6)
b3 <- NormalizeData(b3,verbose = TRUE,scale.factor = 1e6)
VariableFeatures(b1)<-genelist
VariableFeatures(b2)<-genelist
VariableFeatures(b3)<-genelist
datasets <- list(b1,b2,b3)

adata<-merge(datasets[[1]],y=datasets[2:length(datasets)], merge.data=TRUE)
VariableFeatures(adata)<-rownames(adata)
remove(datasets)
adata<-ScaleData(adata)
adata<-RunPCA(adata,npcs = 16)

adata.harmony<-RunHarmony(adata,"X_batch")
saveRDS(adata.harmony, paste0(result_path,"integrated_data.rds"))

