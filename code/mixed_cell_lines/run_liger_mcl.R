library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)

result_path <- "./results/mixed_cell_lines/real/liger/"
dir.create(result_path,showWarnings = FALSE)

infile_b1 <- "./data/mixed_cell_lines/293t_filtered_matrices_mex/hg19/"
infile_b2 <- "./data/mixed_cell_lines/jurkat_filtered_matrices_mex/hg19/"
infile_b3 <- "./data/mixed_cell_lines/mixed_filtered_matrices_mex_permuted/"
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

data.combined <- merge(b1,y=list(b2,b3))
data.combined <- NormalizeData(data.combined,scale.factor = 1e6)
data.combined<-FindVariableFeatures(data.combined, selection.method = "vst", nfeatures = 2000)
# some cell no loading on selected factors of the given genelist.
#VariableFeatures(data.combined) <- genelist
data.combined <- ScaleData(data.combined, split.by = "X_batch", do.center = FALSE)
data.combined <- RunOptimizeALS(data.combined, k = 16, lambda = 5, split.by = "X_batch")
data.combined <- RunQuantileAlignSNF(data.combined, split.by = "X_batch", resolution = 0.3)
data.combined <- RunUMAP(data.combined, dims = 1:ncol(data.combined[["iNMF"]]), reduction = "iNMF")

saveRDS(data.combined, file = paste0(result_path,"liger_integrated_data.rds"))
