
rm(list = ls())
ann_atac_vipcca<-read.csv("./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16/ann_atac.csv",row.names = 1)
# ann_atac_seurat <- readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
ann_atac_vipcca$celltype<-as.character(ann_atac_vipcca$celltype)
ann_atac_vipcca[ann_atac_vipcca$pred_max_score<=0,"celltype"]="unknown"
table(ann_atac_vipcca$celltype)



ann_atac_vipcca_pla <- ann_atac_vipcca[ann_atac_vipcca$celltype=="Platelets",]
ann_atac_vipcca_filtered <- ann_atac_vipcca[ann_atac_vipcca$pred_max_score>0,]
sum(ann_atac_vipcca_pla$pred_max_rawscore>0)
sum(ann_atac_vipcca$pred_max_score>0)

table(ann_atac_vipcca_filtered$celltype)

hist(ann_atac_vipcca$pred_max_score,breaks = 1e5)
g1<-ggplot(ann_atac_vipcca, aes(x=nid,y=pred_max_score))+
	geom_point()
g1
#
rm(list = ls())
celltype.predictions<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
celltype.predictions[celltype.predictions$prediction.score.max<0.5,1]="unknown"
table(celltype.predictions$predicted.id)

