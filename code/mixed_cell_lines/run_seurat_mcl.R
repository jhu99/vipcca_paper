# run Seurat V3 on bipolar
library(Seurat)
library(dplyr)
library(kBET)
library(ggplot2)

rm(list = ls())
result_path <- "./results/mixed_cell_lines/real/seurat/"
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

reference.list <- list(batch1=b1,batch2=b2,batch3=b3)
bipolar.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:16)
bipolar.integrated <- IntegrateData(anchorset = bipolar.anchors, dims = 1:16)

library(ggplot2)
library(cowplot)
DefaultAssay(bipolar.integrated) <- "integrated"
bipolar.integrated <- ScaleData(bipolar.integrated, verbose = FALSE)
bipolar.integrated <- RunPCA(bipolar.integrated, npcs = 16, verbose = FALSE)

bipolar.integrated <- RunUMAP(bipolar.integrated, reduction = "pca", dims = 1:16)

png(paste0(result_path,"seurat_v3_2dplot_batch_umap.png"),width=6, height=4, units="in",res=150)
p1 <- DimPlot(bipolar.integrated, reduction = "umap", group.by = "X_batch")
plot(p1)
dev.off()
png(paste0(result_path,"seurat_v3_2dplot_cell_type_umap.png"),width=6, height=4, units="in",res=150)
p2 <- DimPlot(bipolar.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
							repel = TRUE) + NoLegend()
plot(p2)
dev.off()

saveRDS(bipolar.integrated,file = paste0(result_path,"seurat_v3_integrated_data.rds"))

# write.table(bipolar.integrated@meta.data, file = "./data/bipolar/annotation.txt",
# 						quote=FALSE, col.names = FALSE, sep = "\t")

# bipolar.integrated <- readRDS("./results/bipolar_permutation/seurat/seurat_v3_integrated_data.rds")
# features<-bipolar.integrated@assays$integrated@var.features
# DefaultAssay(bipolar.integrated) <- "RNA"
# unconrrected.data <- subset(bipolar.integrated, features = features)
# unconrrected.data@assays$integrated<-NULL
# saveRDS(unconrrected.data, file = "./results/bipolar_permutation//uncorrected/uncorrected.rds")

# markers <- FindAllMarkers(bipolar.integrated, slot="data", assay = "integrated", 
# 													ident.1="batch1", ident.2="batch2", group.by="")
