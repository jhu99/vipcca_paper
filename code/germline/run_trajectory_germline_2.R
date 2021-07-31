rm(list=ls())
library(cowplot)
library(ggplot2)
library(Seurat)
library(monocle3)
library(rlist)
library(SeuratDisk)
library(patchwork)
library(doParallel)
library("viridis")
setwd("~/github/vipcca_code/")
source("./code/run/evaluation_function.R")
cl <- makeCluster(10)
registerDoParallel(cl)
clusterExport(cl, varlist=c("kbet_for_one_celltype", "get.knn", "kBET"), env=environment())

plot_figures_for_germline<-function(part="female"){
  if(part=="female"){
    #vipcca_res="/VIPCCA/cvae_50_64_16_5/"
    vipcca_res="/VIPCCA/cvae_30_64_16_5/"
    start_cell="F_PGC_4W_embryo2_sc3"
    female_pluripotency_marker<-c("POU5F1","NANOG","KLF4","LIN28A")
    female_PGC_early_marker<-c("KIT","SOX15","DPPA3","ALPL")
    female_Meiosis_marker<-c("SYCP1","SYCE3","PRDM9","SPO11")
    female_PGC_late_marker<-c("DDX4","DAZL","MAEL")
    female_selected_marker<-c("NANOG","POU5F1","DDX4","DAZL")
    markers <- female_selected_marker
    #markers<-c(female_pluripotency_marker,female_PGC_early_marker,female_Meiosis_marker,female_PGC_late_marker)
  }else{
    vipcca_res="/VIPCCA/cvae_20_128_16_3/"
    start_cell="M_PGC_4W_embryo3_sc1"
    male_early_marker<-c("POU5F1","NANOG","KLF4","KIT")
    male_Meiosis_marker<-c("NANOS2","CDK6")
    male_selected_marker<-c("NANOG","POU5F1","CDK6","NANOS2")
    
    markers <- male_selected_marker
    #markers<-c(male_early_marker,male_Meiosis_marker)
  }

  result_dir<-c("uncorrected","desc","harmony","liger","mnn","scanorama","seurat",vipcca_res)
  reduction_names<-c("pca","Embeded_z0.45","harmony","iNMF","pca","scanorama", "pca","vipcca")
  method_names <- c("Before Integration","DESC","Harmony","Liger","MNN","Scanorama", "Seurat V3", "VIPCCA")
  
  g1list<-list()
  g2list<-list()
  mmlist<-c()
  comlist<-c()
  corrall<-c()
  require('FNN')
  #browser()
  for(i in 1:length(result_dir)){
  	result_path=paste0("./results/germline/",part,"/",result_dir[i],"/")
  	print(result_path)
  	datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
  	datasets.integrated <- run_umap_on_seurat(datasets.integrated, reduction=reduction_names[i],min.dist=0.3)
  	#datasets.integrated[,datasets.integrated@meta.data$week=="4W"]@meta.data
  	weeks<-datasets.integrated$week
  	weeks<-as.character(weeks)
  	weeks<-as.numeric(unlist(strsplit(weeks, split="W")))
  	datasets.integrated$weeks<-weeks
  	
  	#DimPlot(datasets.integrated,pt.size=10, cells.highlight = c("F_PGC_4W_embryo2_sc3","M_PGC_4W_embryo3_sc1"))+ggtitle(method_names[i])
  	max.k=300
  	mm <- max.k-MixingMetric(datasets.integrated, grouping.var="batch", reduction=reduction_names[i], dims=1:16, max.k=max.k)
  	mmlist <- c(mmlist,median(mm))
  	datasets.integrated@meta.data["is_all"]="all"
  	# com=run_kbet_for_each_celltype4(datasets.integrated,
  	#                                 reduction = reduction_names[i],
  	#                                 celltype_col = "is_all",
  	#                                 method= method_names[i],
  	#                                 batch_col = "batch")
  	# comlist<-rbind(comlist,com[[1]])
  	# 
  	# 
  	g1<-DimPlot(datasets.integrated,group.by = c("weeks"),pt.size=10)+
  	  scale_color_viridis(discrete = TRUE,option = "D")+ggtitle(method_names[i])
  	# if(i<3)
  	#   g2<-DimPlot(datasets.integrated,group.by = c("batch"),pt.size=10)+ggtitle(method_names[i])
  	# else
    g2<-plot_trajectory_germline(datasets.integrated,result_path,start_cells = start_cell,
  	                             method=method_names[i],reduction = reduction_names[i])
  	g1list<-list.append(g1list,g1)
  	g2list<-list.append(g2list,g2)
  	# 
  	# cell_pseudotime<-pseudotime(readRDS(paste0(result_path,"Li & Guo_cds_pt_on_dr.rds")))
  	# cell_pseudotime<-ifelse(is.infinite(cell_pseudotime),100,cell_pseudotime)
  	# cell_pseudotime<-cell_pseudotime/max(cell_pseudotime)
  	# 
  	# X<-cbind(X,cell_pseudotime[rownames(X)])
  	# corrall<-c(corrall,
  	#            cor.test(cell_pseudotime[colnames(datasets.integrated)],weeks,method = "kendall")$estimate)
  }
  #
  browser()
  #print(corrall)
  # gcorr<-ggplot(data.frame(y=corrall,x=method_names),aes(x,y))+geom_bar(stat = "identity",fill="orange")+
  #   ylab("kendall tau")+theme(axis.title.x = element_blank())
  # # browser()
  # comlist["methods"]=rep(method_names,each=5)
  # colnames(X)<-c(markers,method_names)
  # pseudotime<-c(X[,(length(markers)+1):(length(markers)+8)])
  
  # df_all<-data.frame()
  # for(i in 1:length(markers)){
  #   df<-data.frame(marker=rep(markers[i],dim(X)[1]),exp=rep(X[,i],8),pseudotime=pseudotime,methods=rep(method_names,each=dim(X)[1]))
  #   df_all<-rbind(df_all,df)
  # }
  # # browser()
  # gmarker<-ggplot(data=df_all, aes(x = pseudotime,y = exp,color=methods))+
  #   geom_smooth(se = FALSE,size=5,method = "loess")+theme(text = element_text(size = 60))+ylab("Expression")+
  #   facet_wrap(~marker,nrow = 1)
  # png(paste0("./results/germline/",part,"/marker.png"),width = 20,height = 10,res = 100,units = "in")
  # print(gmarker)
  # dev.off()
  
  g1list<-lapply(g1list, function(x) x+theme_bw(base_size = 60)+theme(axis.title = element_blank()))
  g2list<-lapply(g2list, function(x) x+theme_bw(base_size = 60)+theme(axis.title = element_blank()))
  
  legend_1<-cowplot::get_legend(g1list[[1]]+guides(col=guide_legend("Weeks",nrow = 1))+
  																 	theme(legend.text = element_text(size = 60),
  																 				legend.title = element_text(size = 60),
  																 				legend.direction = "horizontal"))
  legend_2 <-cowplot::get_legend(g2list[[1]]+guides(col=guide_legend("Batch"))+
  																 	theme(legend.text = element_text(size = 60),
  																 				legend.title = element_text(size = 60),
  																 				legend.direction = "horizontal"))
  
  galist<-lapply(g1list, function(x) x+NoLegend())
  gblist<-lapply(g2list, function(x) x+NoLegend())
  
  p1arow<-plot_grid(plotlist=galist,ncol=4)
  p1brow<-plot_grid(plotlist = gblist,ncol=4)
  library(patchwork)
  
  png(paste0("./results/germline/",part,"/traj_by_week.png"),width=60,height =20,unit="in",res=100)
  print((p1arow/legend_1)+plot_layout(ncol = 1,heights = c(5,1)))
  #p1arow
  dev.off()
  png(paste0("./results/germline/",part,"/traj_by_batch.png"),width=60,height =20,unit="in",res=100)
  plot((p1brow/legend_2)+plot_layout(ncol = 1,heights = c(5,1)))
  #p1brow
  dev.off()
  #browser()
  datasets.integrated<-FindNeighbors(datasets.integrated, dims = 1:2,reduction = "umap")
  datasets.integrated<-FindClusters(datasets.integrated,resolution = 0.001)
  DimPlot(datasets.integrated)
  df_mm<-data.frame(method=method_names,mm=mmlist)
  #browser()
  g3<-ggplot(df_mm,aes(x=method,y=mmlist))+
    geom_bar(stat = "identity",fill="orange")+
    ylab("mixing metric")+theme(axis.title.x = element_blank(),text = element_text(size=60))

  g4<-ggplot(comlist,aes(y = ymedian,x =nb))+ylab("KBET")+xlab("% number of samples")+
    geom_line(aes(color=methods),size=10)+theme(text = element_text(size = 60))
  # legend_3<-cowplot::get_legend(g4+theme(legend.text = element_text(size = 60),
  #         legend.title = element_text(size = 60),
  #         legend.direction = "horizontal"))
  png(paste0("./results/germline/",part,"/mixture.png"),width = 20,height = 10,res = 100,units = "in")
  print(g3+g4)
  dev.off()
  #
  
  ## feature plot
  
  result_path=paste0("./results/germline/",part,vipcca_res)
  datasets.integrated<-LoadH5Seurat(paste0(result_path,"output.h5seurat"))
  datasets.integrated <- run_umap_on_seurat(datasets.integrated, reduction="vipcca",min.dist=0.3)
  ## markers
  
  if(part=="female"){
    g5<-FeaturePlot(datasets.integrated,features=female_pluripotency_marker,split.by = "batch",by.col = F,pt.size = 5)
    g6<-FeaturePlot(datasets.integrated,features=female_PGC_early_marker,split.by = "batch",by.col = F,pt.size = 5)
    g7<-FeaturePlot(datasets.integrated,features=female_Meiosis_marker,split.by = "batch",by.col = F,pt.size = 5)
    g8<-FeaturePlot(datasets.integrated,features=female_PGC_late_marker,split.by = "batch",by.col = F,pt.size = 5)
    gc<-g5/g6/g7/g8
  }else{
    g5<-FeaturePlot(datasets.integrated,features=male_early_marker,split.by = "batch",by.col = F,pt.size = 5)
    g6<-FeaturePlot(datasets.integrated,features=male_Meiosis_marker,split.by = "batch",by.col = F,pt.size = 5)
    gc<-g5/g6
  }

  png(paste0("./results/germline/",part,"/featureplot_umap.png"),height=30,width=20,unit="in",res=100)
  print(gc)
  dev.off()
  saveRDS(list(p1arow,p1brow,legend_1,legend_2,g3,g4,gc,gmarker,gcorr),file = paste0("./results/germline/",part,"/ggplots.rds"))
}
debug(plot_figures_for_germline)
plot_figures_for_germline(part = "male")




plot_male<-readRDS(file = "./results/germline/male/ggplots.rds")
plot_figures_for_germline(part = "female")
plot_female <- readRDS(file = "./results/germline/female/ggplots.rds")
plot_main_figure<-function(plot_germ,part="male"){
  g1<-plot_grid(plotlist = plot_germ[c(2,4)],ncol = 1,rel_heights = c(5,1),labels = "a",label_size = 60)
  g2<-plot_grid(plotlist = plot_germ[c(1,3)],ncol = 1,rel_heights = c(5,1),labels = "b",label_size = 60)
  g3<-plot_grid(plotlist = plot_germ[c(5,6)],nrow = 1,rel_widths = c(3,3),labels = c("c","d"),label_size = 60)
  # g4<-plot_germ[[8]]+theme(legend.position = "bottom")
  png(paste0("./results/germline/figures_",part,".png"),width = 60,height = 80, units = "in", res = 100)
  plot(plot_grid(g1,g2,g3,ncol = 1,rel_widths = c(1,1),labels = c("","",""),label_size = 60))
  dev.off()
}
#debug(plot_main_figure)
#undebug(plot_main_figure)
plot_main_figure(plot_male)
plot_main_figure(plot_female,part = "female")



png("./results/germline/figures_corr.png",width = 20,height = 10,units = "in",res = 100)
cor_male<-plot_male[[9]]+ggtitle("Correlation between pseudotime and collection time in male")+
  theme_bw(base_size = 18)+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 45))+labs(tag = "(a)")
cor_female<-plot_female[[9]]+
  ggtitle("Correlation between pseudotime and collection time in female")+
  theme_bw(base_size = 18)+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 45))+labs(tag = "(b)")
cor_male|cor_female
dev.off()
# legend_3<-cowplot::get_legend(plot_male[[6]]+theme(legend.text = element_text(size = 60),
#                                        legend.title = element_text(size = 60),
#                                        legend.direction = "horizontal"))





####



g1<-plot_male[[1]]+draw_label("Male",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
g2<-plot_female[[1]]+draw_label("Female",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
g3<-plot_male[[2]]+draw_label("Male",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
g4<-plot_female[[2]]+draw_label("Female",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
g5<-plot_male[[5]]+ggtitle("Male")+theme(text = element_text(size = 50))+
  plot_female[[5]]+ggtitle("Female")+theme(text = element_text(size = 50))+plot_layout(nrow = 1)
g6<-plot_male[[6]]+ggtitle("Male")+theme(text = element_text(size = 50),legend.position = "none")+ylab("KBET")+xlab("% number of samples")+
  plot_female[[6]]+ggtitle("Female")+theme(text = element_text(size = 50),legend.position = "none")+ylab("KBET")+xlab("% number of samples")+plot_layout(nrow = 1)
legend_methods<-cowplot::get_legend(plot_female[[6]]+theme(legend.title = element_blank(),legend.text = element_text(size=60),legend.direction = "horizontal")+guides(colour = guide_legend(nrow = 1,size=40)))

#plot_male[[5]]
png("./results/germline/traj_by_week.png",width = 60,height = 40, units = "in", res = 100)
(g1/plot_male[[3]]/g2/plot_female[[3]])+plot_layout(ncol = 1,heights = c(4,1,4,1))
dev.off()
png("./results/germline/traj_by_batch.png",width = 60,height = 60, units = "in", res = 100)
g3/plot_male[[4]]/g4/plot_female[[4]]/g5/g6/legend_methods+plot_layout(ncol = 1,heights = c(6,1,6,1,4,4,1))
dev.off()


# df_coef=data.frame(coef=df_coef,method=method_names)
# png("./results/germline/kendall.png",width=20,height=12,unit="in",res=100)
# ggplot(df_coef,aes(x=method,y=coef))+
# 	geom_bar(stat="identity",fill="orange")+theme_bw(base_size=30)+
# 	theme(axis.title.x = element_blank())+ylab("Kendall")
# dev.off()

# #
# rm(list = ls())
# result_path="./results/germline/uncorrected/"
# datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
# dataset.pam<- subset(datasets.integrated, idents = "PAM")
# #dataset.pam<-ScaleData(dataset.pam)
# dataset.lps<- subset(datasets.integrated, idents = "LPS")
# #dataset.lps<-ScaleData(dataset.lps)
# cds.integrated <- readRDS("./results/germline/VIPCCA/cvae_10_128_16/LPS & PAM_cds_pt_on_dr.rds")
# cds.pam <- cds.integrated[,cds.integrated@colData[["condition"]]=="PAM"]
# cds.lps <- cds.integrated[,(cds.integrated@colData[["condition"]]=="LPS")]
# 
# sum(names(pseudotime(cds.pam))==rownames(dataset.pam@meta.data))
# dataset.pam@meta.data$order<-as.factor(rank(pseudotime(cds.pam)))
# dataset.lps@meta.data$order<-as.factor(rank(pseudotime(cds.lps)))
# Idents(dataset.pam)<-"order"
# Idents(dataset.lps)<-"order"
# 
# library(patchwork)
# e1<-DoHeatmap(dataset.pam, group.by = "collection_time",size = 50)+ggtitle("PAM")+NoLegend()+ylab("collection_time")
# e2<-DoHeatmap(dataset.lps, group.by = "collection_time",size = 50)+ggtitle("LPS") +guides(color=FALSE,
#                                                                                           fill = guide_colourbar(barwidth = 5, barheight = 15))+
#   theme(legend.text = element_text(size = 40))
# result_path="./results/germline/"
# result_dir<-c("uncorrected","desc","harmony","mnn","scalign","scanorama","seurat","VIPCCA/s_cvae_10_128_12")
# method_names <- c("Before Integration","DESC","Harmony","MNN","scAlign","Scanorama", "Seurat V3", "VIPCCA")
# ep<-e1+e2
# ep<- ep&theme(title = element_text(size = 60))
# df_coef<-data.frame()
# for(i in 1:length(method_names)){
#   res<-result_dir[i]
#   print(paste0(result_path,res,"/PAM_cds_pt_on_dr.rds"))
#   print(paste0(result_path,res,"/LPS_cds_pt_on_dr.rds"))
#   cds.pam <- readRDS(paste0(result_path,res,"/PAM_cds_pt_on_dr.rds"))
#   cds.lps <- readRDS(paste0(result_path,res,"/LPS_cds_pt_on_dr.rds"))
#   dataset.pam@meta.data$order<-as.factor(rank(pseudotime(cds.pam)))
#   dataset.lps@meta.data$order<-as.factor(rank(pseudotime(cds.lps))) 
#   
#   traj1<-unname(pseudotime(cds.pam))
#   traj2<-unname(pseudotime(cds.lps))
#   t1<-as.integer(substr(pData(cds.pam)[["collection_time"]],1,1))
#   t2<-as.integer(substr(pData(cds.lps)[["collection_time"]],1,1))
#   coef1<-cor.test(traj1,t1,method = "kendall")
#   coef2<-cor.test(traj2,t2,method = "kendall")
#   df_coef<-rbind(df_coef,c(coef1$estimate,coef2$estimate))
#   
#   e3<-DoHeatmap(dataset.pam,label = F,group.by = "order",draw.lines = FALSE)+ylab(method_names[i])+NoLegend()
#   e4<-DoHeatmap(dataset.lps,label = F,group.by ="order", draw.lines = FALSE)+guides(color=FALSE)+NoLegend()
#   ep<-ep+e3+e4
# }
# colnames(df_coef)<-c("pam.coef","lps.coef")
# df_coef$method<-method_names
# ep<-ep&theme(axis.title.y = element_text(size = 50))
# png("./results/germline/heatmap.png",width = 80, height = 60, units = "in", res = 100)
# ep+plot_layout(ncol = 2)
# dev.off()




# 
# # run trajectory for comparison
# #rm(list = ls())
# library(Seurat)
# library(rlist)
# source("./code/run/evaluation_function.R")
# source("./code/germline/cellAlign.R")
# result_path="./results/germline/"
# tlist <- list()
# glist <- list()
# alist <- list()
# result_dir<-c("uncorrected","desc","harmony","mnn","scalign","scanorama","seurat","VIPCCA/cvae_10_128_12")
# reduction_names<-c("pca","Embeded_z0.45","harmony","pca","ALIGNED.CCA","scanorama", "pca","scxx")
# method_names <- c("Before Integration","DESC","Harmony","MNN","scAlign","Scanorama", "Seurat V3", "VIPCCA")
# i=0
# df_corr_global<-data.frame()
# df_corr_local<-data.frame()
# for(res in result_dir){
# 	curr_path=paste0(result_path,res,"/")
# 	print(curr_path)
# 	i=i+1
# 	g2<-plot_global_alignment4(curr_path,method = method_names[i])
# 	g3<-plot_local_alignment(curr_path,method = method_names[i])
# 	glist<-c(glist,g2[1:2],g3[1])
# 	#browser()
# 	df_corr_global <- rbind(df_corr_global, g2[[3]])
# 	df_corr_local <- rbind(df_corr_local, g3[[2]])
# }
# # glist.copy<-glist
# png("./results/germline/corr_global.png",width = 20,height = 20,units = "in", res = 100)
# gexp_global<-ggplot(df_corr_global, aes(x=corr,group=method,color=method))+geom_density(size=3)+theme_bw(base_size = 60)+
# 	theme(legend.position=c(0.4,0.6))+xlab("Correlation")
# gexp_global
# dev.off()
# 
# png("./results/germline/corr_local.png",width = 20,height = 20,units = "in", res = 100)
# gexp_local<-ggplot(df_corr_local, aes(x=corr,group=method,color=method))+geom_density(size=3)+theme_bw(base_size = 60)+
# 	theme(legend.position=c(0.4,0.6))+
# 	xlab("Correlation")
# gexp_local
# dev.off()
# 
# glist<-lapply(glist, function(x) x+theme_bw(base_size = 40))
# legend_stim<-cowplot::get_legend(glist[[2]]+theme(legend.title = element_blank(),
#                                                   legend.text = element_text(size = 60),
#                                                   legend.direction = "horizontal"))
# glist<-lapply(glist, function(x) x+NoLegend())
# 
# library(rlist)
# glist1<-glist[seq(1,24,3)]
# glist2<-glist[seq(2,24,3)]
# glist3<-glist[seq(3,24,3)]
# 
# #glist3<-lapply(glist3, function(x)x+geom_point(size=10))
# 
# for(i in 1:8){
#   glist2[[i]]<-glist2[[i]]+ggtitle(method_names[i])
#   glist3[[i]]<-glist3[[i]]+ggtitle(method_names[i])
# }
# 
# p1<-plot_grid(plotlist=gblist2,ncol = 5)
# p2<-plot_grid(plotlist=glist2,ncol = 5)
# p3<-plot_grid(plotlist = glist3,ncol=5)
# 
# png("./results/germline/all.png",width = 48,height = 45,unit="in",res=100)
# plot_grid(p1,p2,p3,ncol = 1)
# dev.off()
# 
# 
# p1<-p1+ draw_label("PAM",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
# p2<-p2+ draw_label("LPS",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size=60)+theme(plot.margin = margin(80, 0, 0, 0,unit = "pt"))
# p3<-plot_grid(p1,p2,nrow=2,labels = c("a","b"),label_size = 60)
# l1 <- plot_grid(p3,legend_traj,nrow = 2,rel_heights = c(8,1))
# 
# p4<-plot_grid(plotlist = glist3,ncol = 4)
# p4<-p4+ draw_label("Global alignment of trajectories",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
# p4 <- plot_grid(p4,legend_stim,nrow = 2,rel_heights = c(6,1))
# p5<-gexp_global +theme(plot.margin = margin(100, 0, 100, 0,unit = "pt"))
# l2 <- plot_grid(p4,p5,nrow = 1,rel_widths = c(4,2),labels = c("c","d"),label_size = 60)
# 
# p6<-plot_grid(l1,l2,nrow = 2,rel_heights = c(4,2))
# png(paste0(result_path,"traj_inference.png"), width = 6*8, height =6*10, units = "in", res = 100)
# p6
# dev.off()
# 
# mode_pvalues<-compare_correlation(df_coef = df_corr_global)
# mode_pvalues$logpvalue <- -log10(mode_pvalues$pvalue)
# p7<-ggplot(mode_pvalues,aes(x=method,y=logpvalue))+
# 	geom_bar(stat = "identity",fill="orange")+
# 	ylab("-log(p-value)")+theme(axis.title.x = element_blank(),text = element_text(size=40))
# png("./results/germline/corr_gex_pvalue.png",width = 20,height = 20,units = "in",res = 100)
# p7
# dev.off()
