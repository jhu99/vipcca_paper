library(Seurat)
library(SeuratData)
library(SingleCellExperiment)
library(batchelor)
library(Matrix)
# suppressMessages(library(Matrix))
# suppressMessages(library(batchelor))
rm(list = ls())
set.seed(42)
result_path="./results/pancreas/mnn/"
dir.create(result_path, showWarning=FALSE)
pancreas.integrated <- readRDS("./results/pancreas/seurat/seurat_v3_integrated_data.rds")
DefaultAssay(object = pancreas.integrated) <- "RNA"
genes.use <- pancreas.integrated@assays$integrated@var.features
genes.use <- genes.use
pancreas.integrated[["integrated"]] <- NULL

Idents(object = pancreas.integrated) <- "dataset"
celseq <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "celseq")]
celseq2 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "celseq2")]
smartseq2 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "smartseq2")]
fluidigmc1 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "fluidigmc1")]
indrop1 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "indrop1")]
indrop2 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "indrop2")]
indrop3 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "indrop3")]
indrop4 <- pancreas.integrated[genes.use,WhichCells(object = pancreas.integrated, idents = "indrop4")]

mnnc.results <- mnnCorrect(
	as.SingleCellExperiment(celseq),
	as.SingleCellExperiment(celseq2),
	as.SingleCellExperiment(fluidigmc1),
	as.SingleCellExperiment(smartseq2),
	as.SingleCellExperiment(indrop1),
	as.SingleCellExperiment(indrop2),
	as.SingleCellExperiment(indrop3),
	as.SingleCellExperiment(indrop4),
	k=16
)

saveRDS(mnnc.results,"./results/pancreas/mnn/mnnc_results.rds")

mnnc.results<-readRDS("./results/pancreas/mnn/mnnc_results.rds")
mnnc.results<-mnnc.results[,Cells(pancreas.integrated)]
obj<-CreateSeuratObject(assay(mnnc.results),assay = "integrated",meta.data = pancreas.integrated@meta.data)
VariableFeatures(obj)<-rownames(obj)
obj<-ScaleData(object = obj)
obj<-RunPCA(object = obj,npcs = 16,verbose = FALSE)
saveRDS(object = obj, file = "./results/pancreas/mnn/mnn_integrated_data.rds")
