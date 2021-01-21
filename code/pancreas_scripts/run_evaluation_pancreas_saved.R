library(Seurat)
library(cowplot)
library(ggplot2)
library(mclust)
library(cluster)

args = commandArgs(trailingOnly=TRUE)
#result_path=args[1]
result_path="./results/pancreas/"
celltype_col="celltype"
batch_col="dataset"

# standard cvae
datasets.integrated=ReadH5AD(paste0(result_path,"scxx_cvae/output.h5ad"),assay = "integrated")
datasets.integrated@assays$integrated@var.features=datasets.integrated@assays$integrated@data@Dimnames[[1]]
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "scxx", dims = 1:32)
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:32, reduction = "scxx")
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.45)
ari_scxx_cvae <- adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]])
sil_scxx_cvae <- silhouette(as.integer(datasets.integrated[[celltype_col]][[1]]),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
current_path=paste0(result_path,"scxx_cvae/dimplot_")
png(paste0(current_path,"tech_umap.png"),width=8, height=8, units="in",res=150)
p1 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col) + NoLegend() + ggtitle("scxx_cvae")
plot(p1)
dev.off()
png(paste0(current_path,"celltype_umap.png"),width=8, height=8, units="in",res=150)
p2 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE,	repel = TRUE) + NoLegend()+ ggtitle("scxx_cvae")
plot(p2)
dev.off()

# CVAE2, without batch input for encoder
datasets.integrated=ReadH5AD(paste0(result_path,"scxx_cvae2/output.h5ad"),assay = "integrated")
datasets.integrated@assays$integrated@var.features=datasets.integrated@assays$integrated@data@Dimnames[[1]]
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "scxx", dims = 1:32)
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:32, reduction = "scxx")
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.45)
ari_scxx_cvae2 <- adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]])
sil_scxx_cvae2 <- silhouette(as.integer(datasets.integrated[[celltype_col]][[1]]),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
current_path=paste0(result_path,"scxx_cvae2/dimplot_")
png(paste0(current_path,"tech_umap.png"),width=8, height=8, units="in",res=150)
p3 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col)+ NoLegend()+ ggtitle("scxx_cvae2")
plot(p3)
dev.off()
png(paste0(current_path,"celltype_umap.png"),width=8, height=8, units="in",res=150)
p4 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE,	repel = TRUE) + NoLegend()+ ggtitle("scxx_cvae2")
plot(p4)
dev.off()

# scanorama
datasets.integrated=ReadH5AD(paste0(result_path,"scanorama/corrected_subvar.h5ad"),assay = "integrated")
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "scanorama", dims = 1:100)
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:100, reduction = "scanorama")
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.45)
ari_scanorama <- adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]])
sil_scanorama <- silhouette(as.integer(datasets.integrated[[celltype_col]][[1]]),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
png(paste0(result_path,"scanorama/scanorama_dataset_umap.png"),width=8, height=8, units="in",res=150)
p11 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col)+ NoLegend()+ ggtitle("scanorama")
plot(p11)
dev.off()
png(paste0(result_path,"scanorama/scanorama_celltype_umap.png"),width=8, height=8, units="in",res=150)
p12 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE, 
							 repel = TRUE) + NoLegend()+ ggtitle("scanorama")
plot(p12)
dev.off()

# desc
datasets.integrated=ReadH5AD(paste0(result_path,"desc/output.h5ad"),assay = "integrated")
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "Embeded_z0.45", dims = 1:32)
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:32, reduction = "Embeded_z0.45")
ari_desc <- adjustedRandIndex(x=as.integer(datasets.integrated$desc.0.45),y=datasets.integrated[[celltype_col]][[1]])
sil_desc <- silhouette(as.integer(datasets.integrated[[celltype_col]][[1]]),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
png(paste0(result_path,"desc/desc_dataset_umap.png"),width=8, height=8, units="in",res=150)
p13 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col)+ NoLegend()+ ggtitle("desc")
plot(p13)
dev.off()
png(paste0(result_path,"desc/desc_celltype_umap.png"),width=8, height=8, units="in",res=150)
p14 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE, 
							 repel = TRUE) + NoLegend()+ ggtitle("desc")
plot(p14)
dev.off()

# seurat
celltype_col="celltype"
batch_col="dataset"
datasets.integrated<-readRDS(paste0(result_path,"seurat/seurat_v3_integrated_data.rds"))
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:30, reduction = "pca")
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.45)
ari_seurat <- adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]])
sil_seurat <- silhouette(as.numeric(as.factor(datasets.integrated[[celltype_col]][[1]])),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
current_path=paste0(result_path,"seurat/dimplot_")
png(paste0(current_path,"tech_umap.png"),width=8, height=8, units="in",res=150)
p5 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col)+ NoLegend()+ ggtitle("seurat")
plot(p5)
dev.off()
png(paste0(current_path,"celltype_umap.png"),width=8, height=8, units="in",res=150)
p6 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE,	repel = TRUE) + NoLegend()+ ggtitle("seurat")
plot(p6)
dev.off()

# scalign
scAlignobj=readRDS(paste0(result_path,"scalign/scalign_integrated_scalignobj.rds"))
library(dplyr)
scAlignSeuratObj = as.Seurat(scAlignobj, counts="counts", scale.data="scale.data")
datasets.integrated <- RunUMAP(scAlignSeuratObj, reduction = "ALIGNED-MultiCCA", dims = 1:30)
datasets.integrated <- FindNeighbors(datasets.integrated, dims = 1:30, reduction = "ALIGNED-MultiCCA")
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.45)
ari_scalign <- adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]])
sil_scalign <- silhouette(as.numeric(as.factor(datasets.integrated[[celltype_col]][[1]])),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
png(paste0(result_path,"scalign/scalign_dataset_umap.png"),width=8, height=8, units="in",res=150)
p7 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col)+ NoLegend()+ ggtitle("scalign")
plot(p7)
dev.off()
png(paste0(result_path,"scalign/scalign_celltype_umap.png"),width=8, height=8, units="in",res=150)
p8 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE,
							repel = TRUE) + NoLegend()+ ggtitle("scalign")
plot(p8)
dev.off()

# mnn
datasets.integrated <- readRDS(paste0(result_path,"mnn/mnn_integrated_data.rds"))
datasets.integrated <- RunUMAP(datasets.integrated, reduction = "pca", dims = 1:30)
datasets.integrated <- FindClusters(datasets.integrated, resolution = 0.45)
ari_mnn <- adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]])
sil_mnn <- silhouette(as.integer(as.factor(datasets.integrated[[celltype_col]][[1]])),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
png(paste0(result_path,"mnn/mnn_dataset_umap.png"),width=8, height=8, units="in",res=150)
p9 <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col)+ NoLegend()+ ggtitle("mnn")
plot(p9)
dev.off()
png(paste0(result_path,"mnn/mnn_celltype_umap.png"),width=8, height=8, units="in",res=150)
p10 <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col, label = TRUE, 
							 repel = TRUE) + NoLegend()+ ggtitle("mnn")
plot(p10)
dev.off()

# plot grid
png(paste0(result_path,"dataset_umap.png"),width=15, height = 15, units="in", res=100)
plot_grid(p1,p3,p5,p7,p9,p11,p13, labels = c("a","b","c","d","e","f"))
dev.off()
png(paste0(result_path,"celltype_umap.png"),width=15, height = 15, units="in", res=100)
plot_grid(p2,p4,p6,p8,p10,p12,p14, labels = c("a","b","c","d","e","f"))
dev.off()
png(paste0(result_path,"cell_alignment_umap.png"),width=28, height = 8, units="in", res=100)
plot_grid(p1,p3,p5,p7,p9,p11,p13,p2,p4,p6,p8,p10,p12,p14, labels = c("a","b","c","d","e","f","g","h","i","j","k","l"),nrow=2)
dev.off()
# 
# 
png("./results/pancreas/ari.png",width=4, height = 4, units="in", res=100)
barplot(c(ari_scxx_cvae,ari_scxx_cvae2,ari_seurat,ari_scalign,ari_mnn,ari_scanorama,ari_desc),
				ylim = c(0,1.0), names.arg = c("cvae","cvae2","seurat","scalign","mnn","scanorama","desc"),
				ylab = "ARI",las=2)
dev.off()
# Compute or Extract Silhouette Information
data <- data.frame(sil=c(sil_scxx_cvae[,3],sil_scxx_cvae2[,3],
												 sil_seurat[,3],sil_scalign[,3],
												 sil_mnn[,3],sil_scanorama[,3],sil_desc[,3]),
									 group=rep(c("cvae","cvae2","seurat","scalign","mnn","scanorama","desc"),
									 					c(length(sil_scxx_cvae[,3]),length(sil_scxx_cvae[,3]),
									 						length(sil_seurat[,3]),length(sil_scalign[,3]),
									 						length(sil_mnn[,3]),length(sil_scanorama[,3]),
									 						length(sil_desc[,3]))))
png(paste0(result_path,"silhouette.png"),width=4, height = 4, units="in", res=100)
boxplot(sil~group,data=data, ylab="Silhouette Information")
dev.off()

