rm(list=ls())
library(cowplot)
library(ggplot2)
library(Seurat)
library(monocle3)
library(SeuratDisk)
library(rlist)

part="female"
source("./code/run/evaluation_function.R")
# run trajectory on cvae result
result_path=paste0("./results/germline/",part,"/VIPCCA/")
ann_file=paste0("./data/germline/",part,"/annotation.txt")
glist1 <- list()
glist2 <- list()
glist3 <- list()
glist4 <- list()
df_coef<-data.frame()
cvae_list<-c()
for(cvae in dir(result_path, pattern = "cvae")){
  curr_path=paste0(result_path,cvae,"/")
  print(curr_path)
  Convert(paste0(curr_path,"output.h5ad"),dest = "h5Seurat",overwrite = T)
  datasets.integrated<-LoadH5Seurat(paste0(curr_path,"output.h5seurat"))
  datasets.integrated <- run_umap_on_seurat(datasets.integrated, reduction="vipcca",min.dist=0.3)
  g1<-DimPlot(datasets.integrated,group.by = c("batch"),order=c("Guo"))
  g2<-DimPlot(datasets.integrated,group.by = c("batch"),reduction="tsne",order=c("Guo"))
  g1<-g1+ggtitle(cvae)
  glist1 <- list.append(glist1,g1)
  g2<-g2+ggtitle(cvae)
  glist2 <- list.append(glist2,g2)
}

png(paste0(result_path,"all_umap_batch.png"), width = 32, height =32, units = "in", res = 100)
plot_grid(plotlist = glist1, ncol = 5, labels = "auto")
dev.off()
png(paste0(result_path,"all_tsne_batch.png"), width = 32, height =32, units = "in", res = 100)
plot_grid(plotlist = glist2, ncol = 5, labels = "auto")
dev.off()

# feature plot for  result
rm(list=ls())
part="female"
# female FGC and male FGC
source("./code/run/evaluation_function.R")
result_path=paste0("./results/germline/",part,"/VIPCCA/")
# Convert(paste0(result_path,"cvae_50_128_16_5/output.h5ad"),dest = "h5seurat")
ifile=paste0(result_path,"cvae_50_128_16_5/output.h5seurat")
datasets.integrated=LoadH5Seurat(ifile)
datasets.integrated <- run_umap_on_seurat(datasets.integrated, reduction="vipcca",min.dist=0.3)
marker_genes1<-c("POU5F1","NANOG","PRDM9","SYCP1")
marker_genes2<-c("SPO11","DAZL","DDX4","ZP3")
g5<-FeaturePlot(datasets.integrated,features=marker_genes1,split.by = "batch",by.col = F,pt.size = 5)
g6<-FeaturePlot(datasets.integrated,features=marker_genes2,split.by = "batch",by.col = F,pt.size = 5)
png(paste0(result_path,"/featureplot_umap.png"),height=30,width=20,unit="in",res=100)
g5/g6
dev.off()

g8<-FeaturePlot(datasets.integrated,reduction = "tsne",features=marker_genes1,split.by = "batch",by.col = F,pt.size = 5)
g9<-FeaturePlot(datasets.integrated,reduction = "tsne",features=marker_genes2,split.by = "batch",by.col = F,pt.size = 5)
png(paste0(result_path,"/featureplot_tsne.png"),height=30,width=20,unit="in",res=100)
g8/g9
dev.off()

# 
# datasets.integrated <- readRDS("./results/cell_align/uncorrected/integrated_data.rds")
# dataset.pam<- subset(datasets.integrated, idents = "PAM")
# dataset.lps<- subset(datasets.integrated, idents = "LPS")
# medexpress<-c(apply(dataset.pam@assays$integrated@data,2,median),
#               apply(dataset.lps@assays$integrated@data,2,median))
# mdexp1<-apply(dataset.pam@assays$integrated@data,2,median)
# mdexp2<-apply(dataset.lps@assays$integrated@data,2,median)
# 
# curr_path<-"./results/cell_align/uncorrected/"
# cds1<-readRDS(paste0(curr_path,"PAM_cds_pt_on_dr.rds"))
# traj1<-unname(pseudotime(cds1))
# cds2<-readRDS(paste0(curr_path,"LPS_cds_pt_on_dr.rds"))
# traj2<-unname(pseudotime(cds2))
# 
# df<-data.frame(
#   medexpress=medexpress,
#   predtime=c(rescale(traj1),rescale(traj2)),
#   stim=rep(c("PAM","LPS"),times=c(length(traj1),length(traj2)))
# )
# 
# g2<-ggplot(data = df, aes(x=predtime, y=medexpress, colour=stim)) + 
#   geom_point(aes(colour=stim)) + geom_smooth(se = T)+
#   ggtitle("VIPCCA")+ylab("Scaled expression")+xlab("Pseudotime")+
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size=20),
#         legend.position = c(0.4,0.8))
# 
# #######
# rm(list = ls())
# library(cellAlign)
# source("./code/run/evaluation_function.R")
# data(expGlobalLPS)
# data(expGlobalPAM)
# numPts = 200
# curr_path<-"./results/cell_align/VIPCCA/cvae_10_64_16/"
# cds1<-readRDS(paste0(curr_path,"PAM_cds_pt_on_dr.rds"))
# traj1<-rescale(unname(pseudotime(cds1)))
# cds2<-readRDS(paste0(curr_path,"LPS_cds_pt_on_dr.rds"))
# traj2<-rescale(unname(pseudotime(cds2)))
# 
# trajLPS <- pseudotime(cds1)
# trajPAM <- pseudotime(cds2)
# trajLPS<-trajLPS/max(trajLPS)
# trajPAM<-trajPAM/max(trajPAM)
# 
# interGlobalPAM = cellAlign::interWeights(expDataBatch = expGlobalPAM, trajCond = traj1,
#                                          winSz = 0.1, numPts = numPts)
# interGlobalLPS = cellAlign::interWeights(expDataBatch = expGlobalLPS, trajCond = traj2,
#                                          winSz = 0.1, numPts = numPts)
# 
# #scale the interpolated data (Recommended):
# interScaledGlobalLPS = cellAlign::scaleInterpolate(interGlobalLPS)
# interScaledGlobalPAM = cellAlign::scaleInterpolate(interGlobalPAM)
# alignment = globalAlign(interScaledGlobalPAM$scaledData, interScaledGlobalLPS$scaledData,
#                         scores = list(query = interScaledGlobalPAM$traj, 
#                                       ref = interScaledGlobalLPS$traj),
#                         sigCalc = F, numPerm = 20)
# png(paste0(curr_path,"global_alignment_test.png"), width = 800, height = 600)
# plotAlign(alignment)
# dev.off()
# 
# mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalPAM$traj, realTrajQuery = trajPAM,
#                             intTrajRef = interScaledGlobalLPS$traj, realTrajRef = trajLPS)
# 
# png(paste0(result_path,"global_mapping.png"), width = 800, height = 600)
# g2<-plotMapping(mapping)
# print(g2)
# dev.off()
# 
# mdsLPS<-data.frame(mdexp=colMedians(interScaledGlobalLPS$scaledData),
#            traj=(1:200)/200,cond="LPS")
# mdsPAM<-data.frame(mdexp=colMedians(interScaledGlobalPAM$scaledData),
#                    traj=(1:200)/200,cond="PAM")
# 
# mdsALL<-c()
# for(i in 1:length(alignment$align[[1]]$index1)){
#   ppam <- alignment$align[[1]]$index1[i]
#   plps <- alignment$align[[1]]$index2[i]
#   com <- rbind(mdsLPS[plps,],mdsPAM[ppam,])
#   com$paired <- i
#   mdsALL<-rbind(mdsALL,com)
# }
# mdsALL<-mdsALL[c(seq(1,622,2),seq(2,622,2)),]
# s<-seq(1,311,3)
# mdsALL<-mdsALL[c(s,s+311),]
# 
# ggplot(mdsALL,aes(x=traj,y=mdexp,group=cond,color=cond))+
#   geom_point()+geom_line(aes(group=paired),color="grey")+
#   xlab("Scaled pseudotime")+ylab("Scaled expression")
# 
##

# result_path="./results/cell_align2/VIPCCA/"
# glist <- list()
# df_coef<-data.frame()
# cvae_list<-c()
# cvae="cvae_10_128_12"
# for(cvae in dir(result_path, pattern = "cvae")){
# 	curr_path=paste0(result_path,cvae,"/")
# 	print(curr_path)
# 	h5ad_to_seurat(input_file = "output.h5ad", result_path = curr_path,
# 								 annotation_file="./data/cell_align/annotation_dc.txt")
# 	datasets.integrated <- readRDS(paste0(curr_path,"integrated_data.rds"))
# 	g<-plot_trajectory2(datasets.integrated, curr_path)
# 	#browser()
# 
# 	cds <- readRDS(paste0(result_path,cvae,"/LPS & PAM_cds_pt_on_dr.rds"))
# 	cds1 <- cds[,pData(cds)$condition=="LPS"]
# 	cds2 <- cds[,pData(cds)$condition=="PAM"]
# 	traj1<-unname(pseudotime(cds1))
# 	traj2<-unname(pseudotime(cds2))
# 	t1<-as.integer(substr(pData(cds1)[["collection_time"]],1,1))
# 	t2<-as.integer(substr(pData(cds2)[["collection_time"]],1,1))
# 	
# 	coef1<-cor.test(traj1,t1,method = "kendall")
# 	coef2<-cor.test(traj2,t2,method = "kendall")
# 	
# 	df_coef <- rbind(df_coef,c(lps=coef1$estimate,pam=coef2$estimate))
# 	
# 	g<-lapply(g, function(x) return(x+ggtitle(cvae)))
# 	glist <- c(glist,g)
# }
# df_coef$method <- dir(result_path, pattern = "cvae")
# colnames(df_coef)<- c("lps","pam")
# 
# png(paste0(result_path,"all_umap.png"), width = 24*3, height =20*3, units = "in", res = 100)
# plot_grid(plotlist = glist, ncol = 8, labels = "auto")
# dev.off()
