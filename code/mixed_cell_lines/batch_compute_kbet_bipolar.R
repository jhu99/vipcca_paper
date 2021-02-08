library(Seurat)
library(dplyr)
library(kBET)
library(mclust)
library(cluster)

obj <- readRDS("./results/bipolar/scxx_cvae/scxx_cvae_lognorm/scxx_integrated_data.rds")
bipolar.integrated <- obj$data
# kbet alignment score
# png(paste0(result_path,"scxx_kbet_replicate_id.png"),width=6, height=4, units="in",res=150)
# replicate.estimate <- kBET(df = bipolar.integrated@reductions$pca@cell.embeddings,batch = bipolar.integrated@meta.data$replicate.id)
# dev.off()
# png(paste0(result_path,"scxx_kbet_batch.png"),width=6, height=4, units="in",res=150)
# batch.estimate <- kBET(df = bipolar.integrated@reductions$pca@cell.embeddings,batch = bipolar.integrated@meta.data$batch.name)
# dev.off()
# saveRDS(list(replicate.estimate=replicate.estimate,batch.estimate=batch.estimate),file = paste0(result_path,"scxx_kbet_score.rds"))
# Clustering quality
sil_scxx <- silhouette(as.integer(bipolar.integrated$seurat_clusters),dist = dist(bipolar.integrated@reductions$pca@cell.embeddings))
sil_scxx.summ <- summary(sil_scxx)
ari_scxx <- adjustedRandIndex(x=as.integer(bipolar.integrated$seurat_clusters),y=bipolar.integrated$cell.type)


obj <- readRDS("./results/bipolar/seurat/")
bipolar.integrated <- obj$data
# kbet alignment score
# png(paste0(result_path,"seurat_v3_kbet_replicate_id.png"),width=6, height=4, units="in",res=150)
# replicate.estimate <- kBET(df = bipolar.integrated@reductions$pca@cell.embeddings,batch = bipolar.integrated@meta.data$replicate_id)
# dev.off()
# png(paste0(result_path,"seurat_v3_kbet_batch.png"),width=6, height=4, units="in",res=150)
# batch.estimate <- kBET(df = bipolar.integrated@reductions$pca@cell.embeddings,batch = bipolar.integrated@meta.data$batch_name)
# dev.off()
# saveRDS(list(replicate.estimate=replicate.estimate,batch.estimate=batch.estimate),file = paste0(result_path,"seurat_v3_kbet_score.rds"))
sil_seurat <- silhouette(as.integer(bipolar.integrated$seurat_clusters),dist = dist(bipolar.integrated@reductions$pca@cell.embeddings))
sil_seurat.summ <- summary(sil_seurat)

# plot silhouette
df <- list(clus.avg.widths=c(sil_scxx.summ$clus.avg.widths, sil_seurat.summ$clus.avg.widths), method=rep(c("scxx","seurat"), each=15))
png(paste0(result_path,"evaluation_cluster_silhouette.png"),width=6, height=4, units="in",res=150)
boxplot(clus.avg.widths~method, data=df,ylab="silhouette")
dev.off()

ari_seurat <- adjustedRandIndex(x=as.integer(bipolar.integrated$seurat_clusters),y=bipolar.integrated$celltype)
png(paste0(result_path,"evaluation_cluster_ari.png"),width=6, height=4, units="in",res=150)
height=c(ari_scxx,ari_seurat)
mp<-barplot(height = height,names.arg = c("scxx","seurate V3"),ylim = 0:1,ylab = "ARI")
text(mp, height, labels = format(height, 4), pos = 1)
dev.off()
