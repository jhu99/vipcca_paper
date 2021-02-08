library(cowplot)
library(ggplot2)
library(doParallel)
library(kBET)
library(FNN)
library(gridExtra)
library(grid)
library(Seurat)

rm(list = ls())
set.seed(2020)

source("./code/run/evaluation_function.R")
result_path="./results/mixed_cell_lines/permuted/"
celltype_col="celltype"
batch_col="X_batch"
#cvae
curr_path=paste0(result_path,"scxx_cvae/cvae_2_64_8/")
cvae_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="VIPCCA",algorithm="cover_tree",
																			 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
curr_path=paste0(result_path,"liger/")
liger_res <- calculate_alignment_metric(curr_path,
																				"liger_integrated_data.rds",method="liger",
																				batch_col=batch_col,celltype_col=celltype_col,reduction="iNMF",dims=1:20,fmt="rds")

curr_path=paste0(result_path,"seurat/")
seurat_res <- calculate_alignment_metric(curr_path,"seurat_v3_integrated_data.rds",method="seurat",
																				 batch_col=batch_col,celltype_col=celltype_col,reduction="pca",dims=1:16,fmt="rds")

curr_path=paste0(result_path,"scanorama/")
scanorama_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="scanorama",
																						batch_col=batch_col,celltype_col=celltype_col,reduction="scanorama",dims=1:16,fmt="rds")
library(rlist)
ga<-list()
gb<-list()
gd<-list()
ari<-c()
sil<-c()
mm<-c()
tt<-c()
rr<-c()
lrr<-c()
for(res in list(cvae_res,liger_res,seurat_res,scanorama_res)){
	ga<-list.append(ga,res[[4]])
	#
	mm <- c(mm,res[[1]])
	sil <- c(sil,res[[2]])
	ari <- c(ari,res[[3]])
	tt<- c(tt,length(res[[2]]))
	rr<-c(rr,res[[6]][[1]][,3])
	lrr<-c(lrr,length(res[[6]][[1]][,3]))
}
for(res in list(cvae_res,liger_res,seurat_res,scanorama_res)){
	ga<-list.append(ga,res[[5]])
	gd<-list.append(gd,res[[6]][2])
}
# ari_data<- data.frame(cari=ari,group=rep(c("cvae",
# 																					 "liger","seurat","scanorama"),
# 																				 each=10))
# 
# sil_data<- data.frame(csil=sil,group=rep(c("cvae",
# 																					 "liger","seurat","scanorama"),
# 																				 tt))
# kbet_data<- data.frame(mkb=rr,group=rep(c("cvae",
# 																					"liger","seurat","scanorama"),																				lrr),
# 											 neighbor_size=rep(seq(0.15,0.95,by = 0.05),2*4),
# 											 cell_type=rep(rep(c("293t","Jurkat"),each=17),4))
# 
# # plot metric
# g1<-ggplot(data = ari_data, aes(x=group,y=cari))+
# 	geom_boxplot()+
# 	ylab("ARI")+
# 	theme(axis.title.x = element_blank())
# 
# # Compute or Extract Silhouette Information
# g2<-ggplot(data = sil_data, aes(x=group,y=csil))+
# 	geom_boxplot()+
# 	ylab("Silhouette")+
# 	theme(axis.title.x = element_blank())
# 
# df=data.frame(mm=mm,aa=c("cvae",
# 												 "liger","seurat","scanorama"))
# g3<-ggplot(data = df,aes(x=aa,y=mm))+
# 	geom_bar(stat="identity")+
# 	ylab("Mixing Metric")+
# 	theme(axis.title.x = element_blank())
# 
# g4<-ggplot(data = kbet_data, aes(x=group,y=mkb))+
# 	geom_boxplot()+
# 	ylab("kBET (acceptance rate)")+
# 	theme(axis.title.x = element_blank())
# 
# # png(paste0(result_path,"kbet_celltype.png"),width=8, height = 16, units="in", res=100)
# # plot_grid(plotlist= gd, labels = "auto", ncol = 1)
# # dev.off()
# kbet_data_sub = kbet_data[kbet_data$neighbor_size==0.2,]
# g5<-ggplot(data= kbet_data_sub,aes(x=cell_type,y = mkb,fill=group))+
# 	geom_bar(position="dodge",stat="identity")+
# 	ylab("kBET acceptance rate")+
# 	xlab("cell types")+
# 	theme(legend.title = element_blank(),
# 				legend.position = c(0.9,0.9))
# # legend3 <- cowplot::get_legend(g5)
# # g5<-g5+theme(legend.position = "none")
# 
# curr_dirs<-c(
# 	"scanorama/integrated_data.rds",
# 	"seurat/seurat_v3_integrated_data.rds",
# 	# "liger/liger_integrated_data.rds",
# 	"scxx_cvae/cvae_2_64_8/integrated_data.rds"
# )
# br.all<-data.frame()
# for(cdir in curr_dirs){
# 	rfile<-paste0(result_path,cdir)
# 	if(file.exists(rfile))
# 		print(rfile)
# 	data.integrated <- readRDS(rfile)
# 	# for each cell type
# 	# foreach(ct=levels(factor(data.integrated@meta.data[[celltype_col]])),
# 	# 				.export = 'run_dge_analyse',
# 	# 				.combine = 'rbind')%dopar% run_dge_analyse(data.integrated,
# 	# 																									 batch1 = "indrop3",
# 	# 																									 batch2 = "smartseq2",
# 	# 																									 celltype1 = ct,
# 	# 																									 celltype_col = celltype_col,
# 	# 																									 batch_col = batch_col,
# 	# 																									 method = strsplit(cdir,split = "/")[[1]][1])
# 	# 
# 	for(ct in levels(factor(data.integrated@meta.data[[celltype_col]]))
# 	){
# 		br<-run_dge_analyse(data.integrated,
# 												batch1 = "mixed",
# 												batch2 = ct,
# 												celltype1 = ct,
# 												celltype_col = celltype_col,batch_col = batch_col,method = strsplit(cdir,split = "/")[[1]][1])
# 		if(!is.null(br)){
# 			br.all<-rbind(br.all,br)
# 		}
# 		else{
# 			print("Cells are less than 50 in at least one group.")
# 		}
# 		print(ct)
# 	}
# }
# saveRDS(br.all,"./results/mixed_cell_lines/real/result_R/sig_genes_all.rds")
# br.all<-readRDS("./results/mixed_cell_lines/real/result_R/sig_genes_all.rds")
# br.all<-br.all[br.all$method!="liger",]
# br.all<-br.all[br.all$method!="desc",]
# br.selected <- br.all[br.all$p_val_adj<1e-20,]
# sig_genes<-as.data.frame(table(br.selected$method,br.selected$ctype))
# 
# g6<-ggplot(data = sig_genes,aes(x = Var1, y=Freq,fill=Var2))+
# 	geom_bar(position="stack", stat="identity")+
# 	ylab("# significant genes")+theme(axis.title.x = element_blank())+
# 	labs(fill = "cell type")+theme(legend.position = c(0.1,0.9))

gx<-list.append(ga)
png(paste0(result_path,"plot.all.png"),width = 20, height = 10, units = "in", res=200)
plot_grid(plotlist = gx,ncol = 4,labels="auto")
dev.off()
