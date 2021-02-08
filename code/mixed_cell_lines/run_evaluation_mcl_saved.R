library(cowplot)
library(ggplot2)
library(doParallel)
library(kBET)
library(FNN)
library(gridExtra)
library(grid)
library(Seurat)

set.seed(2020)
rm(list = ls())
source("./code/run/evaluation_function.R")
result_path="./results/mixed_cell_lines/real/"
celltype_col="celltype"
batch_col="X_batch"
dims=1:16
min.dist=0.01
reduction="scxx"


dims=1:16
reduction="scxx"
resolution=0.2
h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(result_path,"scxx_cvae/cvae_1_64_8/"),
							 annotation_file="./data/mixed_cell_lines/annotation.txt",
#							 col.names=c("barcode","orig.ident","nCount_RNA","nFeature_RNA","celltype"
#							 						,"cluster","replicate_id","batch_name")
)
h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(result_path,"scanorama/"),
							 annotation_file="./data/mixed_cell_lines/annotation_scanorama.txt",
							 #							 col.names=c("barcode","orig.ident","nCount_RNA","nFeature_RNA","celltype"
							 #							 						,"cluster","replicate_id","batch_name")
)
curr_path=paste0(result_path,"scxx_cvae/cvae_2_64_8/")
datasets.cvae<-readRDS(paste0(curr_path,"integrated_data.rds"))

cvae_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="VIPCCA",algorithm="cover_tree",
																			 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")


datasets.cvae <- RunUMAP(datasets.cvae, reduction = reduction, dims = dims, min.dist = 0.1, n_threads=10, scale=TRUE)
datasets.cvae <- FindNeighbors(datasets.cvae, dims = dims, reduction = reduction)
datasets.cvae <- FindClusters(datasets.cvae, resolution = resolution)
g1<-DimPlot(datasets.cvae,reduction = "umap", group.by = batch_col)+ NoLegend()+ggtitle("VIPCCA")
g2<-DimPlot(datasets.cvae,reduction = "umap", group.by = celltype_col)+NoLegend()+ggtitle("VIPCCA")
c1<-DimPlot(datasets.cvae,reduction = "umap")+ggtitle("cvae")
method="cvae"
png(paste0(curr_path,"plot_all.png"),width = 1200, height = 400)
plot_grid(g1,g2,c1,nrow = 1)
dev.off()



com.cvae=run_kbet_for_each_celltype2(datasets.integrated = datasets.cvae,
																reduction = reduction,
																celltype_col = "celltype",
																batch_col = batch_col,
																method = method,
																algorithm="cover_tree")
# saveRDS(com.cvae,file = "./results/mixed_cell_lines/scxx_cvae/scmb_cvae_1_64_8/kbet.rds")
# remove(datasets.cvae)

reduction="pca"
dims=1:16
resolution=0.5
datasets.seurat<-readRDS(paste0(result_path,"seurat/seurat_v3_integrated_data.rds"))
datasets.seurat <- FindNeighbors(datasets.seurat, dims = dims, reduction = reduction)
datasets.seurat <- FindClusters(datasets.seurat, resolution = resolution)
g3<-DimPlot(datasets.seurat,reduction = "umap", group.by = batch_col)
g4<-DimPlot(datasets.seurat,reduction = "umap", group.by = celltype_col)
legend1 <- cowplot::get_legend(g3)
legend2 <- cowplot::get_legend(g4)
g3<-g3+NoLegend()+ggtitle("seurat")
g4<-g4+NoLegend()+ggtitle("seurat")
c2<-DimPlot(datasets.seurat,reduction = "umap")+ggtitle("seurat")
method="seurat"
reduction="pca"
png(paste0(result_path,"seurat/plot_all.png"), width = 1200, height = 400)
plot_grid(g3,g4,c2, nrow = 1)
dev.off()
# com.seurat=run_kbet_for_each_celltype2(datasets.integrated = datasets.seurat,
# 																reduction = reduction,
# 																celltype_col = "seurat_clusters",
# 																batch_col = batch_col,
# 																method = method,
# 																algorithm="cover_tree")
# saveRDS(com.seurat,paste0(result_path,"seurat/kbet.rds"))
# remove(datasets.seurat)

method="liger"
reduction="iNMF"
resolution=0.5
dims=1:16
datasets.liger<-readRDS(paste0(result_path,"/liger/liger_integrated_data.rds"))
g5<-DimPlot(datasets.liger,reduction = "umap", group.by = batch_col)+NoLegend()+ggtitle("liger")
g6<-DimPlot(datasets.liger,reduction = "umap", group.by = celltype_col)+NoLegend()+ggtitle("liger")
datasets.liger <- FindNeighbors(datasets.liger, dims = dims, reduction = reduction)
datasets.liger <- FindClusters(datasets.liger, resolution = resolution)
c3 <- DimPlot(datasets.liger, reduction="umap")+ggtitle("liger")
png(paste0(result_path,"liger/plot_all.png"), width = 1200, height = 400)
plot_grid(g5,g6,c3, nrow = 1)
dev.off()


# com.liger=run_kbet_for_each_celltype2(datasets.integrated = datasets.liger,
# 																reduction = reduction,
# 																celltype_col = "seurat_clusters",
# 																batch_col = batch_col,
# 																method = method,
# 																algorithm="cover_tree")
# remove(datasets.liger)
# saveRDS(com.liger,file = paste0(result_path,"liger/kbet.rds"))

dims=1:16
reduction="scanorama"
resolution=0.5
batch_col="X_batch"
celltype_col="celltype"
datasets.scanorama<- readRDS(paste0(result_path,"scanorama/integrated_data.rds"))
datasets.scanorama <- RunUMAP(datasets.scanorama, reduction = reduction, dims = dims, min.dist = 0.01, n_threads=10, scale=TRUE)
datasets.scanorama <- FindNeighbors(datasets.scanorama, dims = dims, reduction = reduction)
datasets.scanorama <- FindClusters(datasets.scanorama, resolution = resolution)
g7<-DimPlot(datasets.scanorama,reduction = "umap", group.by = batch_col)+NoLegend()+ggtitle("scanorama")
g8<-DimPlot(datasets.scanorama,reduction = "umap", group.by = celltype_col)+NoLegend()+ggtitle("scanorama")
c4<-DimPlot(datasets.scanorama,reduction = "umap")+ggtitle("scanorama")
method="scanorama"
png(paste0(result_path,"scanorama/plot_all.png"),width = 1200, height = 400)
plot_grid(g7,g8,c4,nrow = 1)
dev.off()




g<-plot_grid(plotlist = list(g1,g3,g5,g7,g2,g4,g6,g8),nrow = 2)
p<-plot_grid(plotlist= list(g,plot_grid(legend1,legend2,ncol = 1)), rel_widths = c(3,1))
x<-plot_grid(plotlist = list(c1,c2,c3,c4),  nrow = 1)

png(paste0(result_path,"plot_all.png"),width = 20,height = 8,units = "in",res=100)
print(p)
dev.off()

png(paste0(result_path,"plot_cluster.png"), width = 16, height = 4, units = "in", res=100)
print(x)
dev.off()

#z<-plot_grid(plotlist=list(com.cvae[[2]],com.seurat[[2]],com.liger[[2]]),nrow = 3)

# png("./results/bipolar_permutation/plot_kbet.png",width = 45, height = 9, units = "in", res = 100)
# print(z)
# dev.off()

# kb<-rbind(com.cvae[[1]],com.seurat[[1]],com.liger[[1]])
# kb$method <- rep(c("cvae","seurat","liger"),each=255)
# 
# yx <- ggplot(data=kb,aes(x=method,y=ymedian))+geom_boxplot()
# png("././results/bipolar_permutation/plot_kbet_boxplot.png",width = 6, height = 4, units = "in", res=100)
# plot(yx)
# dev.off()

