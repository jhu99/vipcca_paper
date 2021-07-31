rm(list=ls())
library(cowplot)
library(ggplot2)
library(Seurat)
library(rlist)
library(patchwork)
library(SingleCellExperiment)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library("viridis")
source("./code/run/evaluation_function.R")
part="male"
vipcca<-ifelse(part=="male","/VIPCCA/cvae_20_128_16_3/","/VIPCCA/cvae_30_64_16_5/")
respath<-ifelse(part=="male","./results/germline/male/","./results/germline/female/")
lpos<-ifelse(part=="male",0.3,0.5)

result_dir<-c("uncorrected","desc","harmony","liger","mnn","scanorama","seurat",vipcca)
reduction_names<-c("pca","Embeded_z0.45","harmony","iNMF","pca","scanorama", "pca","vipcca")
method_names <- c("Before Integration","DESC","Harmony","Liger","MNN","Scanorama", "Seurat V3", "VIPCCA")
g1list<-list()
g2list<-list()
for(i in 1:length(result_dir)){
  filename<-paste0(respath,result_dir[i],"/integrated_data.rds")
  datasets.integrated <- readRDS(filename)
  datasets.integrated <- run_umap_on_seurat(datasets.integrated, reduction=reduction_names[i],min.dist=0.3)
  weeks<-datasets.integrated$week
  weeks<-as.character(weeks)
  weeks<-as.numeric(unlist(strsplit(weeks, split="W")))
  datasets.integrated$weeks<-weeks
  g1<-DimPlot(datasets.integrated,group.by = c("weeks"),pt.size=10)+
    scale_color_viridis(discrete = TRUE,option = "D")+ggtitle(method_names[i])
  g2 <- DimPlot(datasets.integrated,group.by = c("batch"),pt.size=10)+ggtitle(method_names[i])
  # slingshot
  sce<-as.SingleCellExperiment(datasets.integrated)
  set.seed(1)
  rd1<-reducedDim(sce,"UMAP")
  colData(sce)$kmeans <- kmeans(rd1,centers = 3,iter.max=100000,nstart=5000)$cluster
  sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'UMAP')
  x<-as.data.frame(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ])
  g1<-g1+geom_point(data = x,mapping = aes(x=UMAP_1, y=UMAP_2,size=3),show.legend = F)
  g1list <- list.append(g1list,g1)
  g2list <- list.append(g2list,g2)
  # mm
  sce$slingPseudotime_1
}
galist<-lapply(g1list, function(x) x+NoLegend()+
                 theme(title = element_text(face = "bold",size=70),axis.title = element_blank(),axis.text  = element_blank(),axis.ticks = element_blank()))
gblist<-lapply(g2list, function(x) x+NoLegend()+
                 theme(title = element_text(face = "bold",size=70),
                       axis.title = element_blank(),axis.text  = element_blank(),axis.ticks = element_blank()))
legend_1<-cowplot::get_legend(g1list[[1]]+
                                theme(legend.text = element_text(size = 70,face = "bold"),
                                      legend.title = element_text(size = 70,face="bold"),
                                      legend.direction = "horizontal",legend.position = c(lpos,0.5))+
                                guides(col=guide_legend("Weeks",nrow = 1,override.aes = list(size=20))))
legend_2 <-cowplot::get_legend(g2list[[1]]+
                                 theme(legend.text = element_text(size = 70,face="bold"),
                                       legend.title = element_text(size = 70,face="bold"),
                                       legend.direction = "horizontal",legend.position = c(0.4,0.5))+
                                 guides(col=guide_legend("Batch",override.aes = list(size=20))))
row1<-(plot_grid(plotlist=galist,ncol=4)/legend_1)+plot_layout(ncol = 1,heights = c(5,1))
row2<-(plot_grid(plotlist=gblist,ncol=4)/legend_2)+plot_layout(ncol = 1,heights = c(5,1))
png(paste0(respath,"traj_by_week.png"),width=40,height =20,unit="in",res=100)
print(row1)
dev.off()
png(paste0(respath,"traj_by_batch.png"),width=40,height =20,unit="in",res=100)
print(row2)
dev.off()
#
mmlist<-list()
comlist<-list()
for(i in 1:length(result_dir)){
  filename<-paste0(respath,result_dir[i],"/integrated_data.rds")
  datasets.integrated <- readRDS(filename)
  max.k=300
  mm <- max.k-MixingMetric(datasets.integrated, grouping.var="batch", 
                           reduction=reduction_names[i], dims=1:16, max.k=max.k)
  mmlist <- c(mmlist,median(mm))
  datasets.integrated@meta.data["is_all"]="all"
  com=run_kbet_for_each_celltype4(datasets.integrated,
                                  reduction = reduction_names[i],
                                  celltype_col = "is_all",
                                  method= method_names[i],
                                  batch_col = "batch")
  comlist<-rbind(comlist,com[[1]])
}
# colors=c("cyan","chocolate","coral","chartreuse",
#          "burlywood","blue","orchid","deeppink1")
# ScAlign out, Before Integration in
colors=c("grey","cyan","chocolate","coral","chartreuse","blue","orchid","deeppink1")
df_mm<-data.frame(method=method_names,mm=mmlist)
g3<-ggplot(df_mm,aes(x=method,y=mmlist))+
  geom_bar(stat = "identity",fill=colors)+ylab("Mixing metric")+
  theme(axis.title.x = element_blank(),text = element_text(size=70,face = "bold"),
        panel.background = element_blank(),axis.text.x = element_blank())
#,axis.text.x = element_text(angle = 45,vjust = 0.5)
comlist["methods"]=rep(method_names,each=5)
g4<-ggplot(comlist,aes(y = ymedian,x =nb))+ylab("KBET")+xlab("% Sample size")+
  geom_line(aes(color=methods),size=10)+
  scale_color_manual(values=colors)+
  theme(text = element_text(size = 70,face="bold"),panel.background = element_blank(),legend.text = element_text(size=60),
        legend.title = element_blank(),legend.direction = "horizontal",legend.position = "bottom",
        legend.key.size=unit(4,"cm"),legend.key.width=unit(4,"cm"))+guides(col=guide_legend(nrow = 2))
legend1<-cowplot::get_legend(g4)
g4<-g4+NoLegend()
row3<-plot_grid(g3,g4,nrow = 1,labels = c("c","d"),label_size = 100,label_x = 0, label_y = 1,
                hjust = 0.5, vjust = 0.5)+theme(plot.margin = unit(c(1,0,0,1),units = "cm"))
png(paste0(respath,"row3.png"),width=50,height =20,unit="in",res=100)
print(row3)
dev.off()
png(paste0(respath,"all-in-one.png"),width = 50,height = 70,units = "in",res = 100)
plot_grid(plotlist = list(row2,row1,row3,legend1),ncol = 1,labels = c("a","b"),label_size = 100,label_x = 0, label_y = 1,
          hjust = 0.5, vjust = 0.5,rel_heights = c(2,2,2,0.5))+theme(plot.margin = unit(c(1,0,0,1),units = "cm"))
dev.off()

##

#sce$slingPseudotime_1




##
# rm(list=ls())
# library(SingleCellExperiment)
# library(slingshot)
# library(grDevices)
# library(RColorBrewer)
# source("./code/run/evaluation_function.R")
# 
# datasets.integrated <- readRDS("./results/germline/male/harmony/integrated_data.rds")
# datasets.integrated <- run_umap_on_seurat(datasets.integrated, reduction="harmony",min.dist=0.3)
# sce<-as.SingleCellExperiment(datasets.integrated)
# set.seed(1)
# rd1<-reducedDim(sce,"UMAP")
# colData(sce)$kmeans <- kmeans(rd1,centers = 3,iter.max=100000,nstart=5000)$cluster
# sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'UMAP')
# x<-as.data.frame(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ])
# g3<-DimPlot(datasets.integrated,group.by = c("weeks"),pt.size=10)+
#   geom_point(data = x,mapping = aes(x=UMAP_1, y=UMAP_2))



# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
# plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sce),lwd=2, col='black')
# 
# colnames(x)
# ggplot(x, aes(x=UMAP_1, y=UMAP_2)) + geom_point()
# 


# rdslits<-dir("./results/germline/",pattern = "cds_expressdata.rds",recursive = T)
# for(rds in rdslits){
#   cmd<-paste0("sshpass -p jhu332 scp ", rds, " jhu@10.69.8.219:~/github/vipcca_code/results/germline/",rds)
#   print(cmd)
# }
