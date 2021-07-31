library(cowplot)
library(ggplot2)
library(doParallel)
library(kBET)
library(FNN)
library(gridExtra)
library(grid)

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
source("./code/run/evaluation_function.R")
result_path="./results/pancreas/"
batch_col="dataset"
celltype_col="celltype"
# browser()
# lambda_regulizer=args[1]
# batch_input_size=args[2]
cl <- makeCluster(20)
registerDoParallel(cl)
clusterExport(cl, varlist=c("kbet_for_one_celltype", "get.knn", "kBET"), env=environment())
# res_list <- list()
# ari<-c()
# gari<-c()
# rr<-c()
# lrr<-c()
# gd<-list()
# curr_path=paste0(result_path,"scxx_cvae/")
# if(!file.exists(paste0(curr_path,"kbet_cvae.png"))){
# 	methods<-dir(curr_path,pattern = "var2_5")
# 	for(dir_cvae in methods){
# 		print(dir_cvae)
# 		h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(curr_path,dir_cvae,"/"),
# 									 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 									 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# 		cvae_res <- calculate_alignment_metric(paste0(curr_path,dir_cvae,"/"),"integrated_data.rds",method=dir_cvae,algorithm="cover_tree",
# 																								batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
# 		ari <- rbind(ari,cvae_res[[3]])
# 		gari <- c(gari,length(cvae_res[[3]]))
# 		rr<-c(rr,cvae_res[[6]][[1]][,3])
# 		lrr<-c(lrr,length(cvae_res[[6]][[1]][,3]))
# 		gd[[dir_cvae]] <- cvae_res[[6]][[2]]
# 	}
# 	ari_data_cvae<- data.frame(cari=ari,group=rep(methods, each=10))
# 	kbet_data_cvae<- data.frame(mkb=rr,group=rep(methods,lrr))
# 	png(paste0(curr_path,"kbet_cvae.png"),width=12, height = 12, units="in", res=100)
# 	boxplot(mkb~group, data=kbet_data_cvae, ylab="kBET (acceptance rate)", las=2)
# 	dev.off()
# 	png(paste0(curr_path,"ari_cvae.png"),width=12, height = 12, units="in", res=100)
# 	boxplot(cari~group,data = ari_data_cvae, ylab="Adjusted Rand Index", las=2)
# 	dev.off()
# }
# 
# saveRDS(list(ari_data_cvae,kbet_data_cvae),file = "./results/pancreas/scxx_cvae_1e4/cvae_stat.rds")
# 
# # head(kbet_data_cvae)
# # head(ari_data_cvae)
# # saveRDS(list(kbet_data_cvae,ari_data_cvae),file = paste0(curr_path,"cvae_stat_1235.rds"))
# #rm(list = ls())
# res_list <- list()
# ari<-c()
# gari<-c()
# rr<-c()
# lrr<-c()
# gd<-list()
# curr_path=paste0(result_path,"scxx_cvae2/")
# if(!file.exists(paste0(curr_path,"kbet_cvae2.png"))){
# 	methods<-dir(curr_path,pattern = "cvae2_2")
# 	for(dir_cvae in methods){
# 		print(dir_cvae)
# 		h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(curr_path,dir_cvae,"/"),
# 									 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 									 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# 		cvae_res <- calculate_alignment_metric(paste0(curr_path,dir_cvae,"/"),"integrated_data.rds",method=dir_cvae,algorithm="cover_tree",
# 																					 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
# 		ari <- c(ari,cvae_res[[3]])
# 		gari <- c(gari,length(cvae_res[[3]]))
# 		rr<-c(rr,cvae_res[[6]][[1]][,3])
# 		lrr<-c(lrr,length(cvae_res[[6]][[1]][,3]))
# 		gd[[dir_cvae]] <- cvae_res[[6]][[2]]
# 	}
# 	ari_data_cvae2<- data.frame(cari=ari,group=rep(methods, each=10))
# 	kbet_data_cvae2<- data.frame(mkb=rr,group=rep(methods,lrr))
# 	png(paste0(curr_path,"kbet_cvae.png"),width=12, height = 12, units="in", res=100)
# 	boxplot(mkb~group, data=kbet_data_cvae2, ylab="kBET (acceptance rate)", las=2)
# 	dev.off()
# 	png(paste0(curr_path,"ari_cvae.png"),width=12, height = 12, units="in", res=100)
# 	boxplot(cari~group,data = ari_data_cvae2, ylab="Adjusted Rand Index", las=2)
# 	dev.off()
# }
# saveRDS(list(ari_data_cvae2,kbet_data_cvae2),file = "./results/pancreas/scxx_cvae2/cvae2_stat_1510.rds")
# ax1<-readRDS("./results/pancreas/scxx_cvae/cvae_stat.rds")
# ax2<-readRDS("./results/pancreas/scxx_cvae2/cvae2_stat_1510.rds")
# ax3<-readRDS("./results/pancreas/scxx_cvae_1e4/cvae_stat.rds")
# ad1<-aggregate(cari~group,ax1[[1]],max)
# md1<-aggregate(mkb~group,ax1[[2]],median)
# md1[md1$mkb>0.9,]
# ad1[rownames(md1[md1$mkb>0.9,]),]
# 
# ad2<-aggregate(cari~group,ax2[[1]],max)
# md2<-aggregate(mkb~group,ax2[[2]],median)
# md2[md2$mkb>0.8,]
# ad2[rownames(md2[md2$mkb>0.8,]),]
# 
# ad3<-aggregate(cari~group,ax3[[1]],mean)
# md3<-aggregate(mkb~group,ax3[[2]],median)
# md3[md3$mkb>0.9,]
# ad3[rownames(md3[md3$mkb>0.9,]),]


#saveRDS(res, "./results/pancreas/scxx_cvae_res.rds")
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scxx_cvae/scxx_cvae_l1_b128/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scxx_cvae2/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scanorama/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/desc/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scvi/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scgen/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
vipcca_res <- calculate_alignment_metric(paste0(result_path, "scxx_cvae/cvae_var2_10_128_16/"),"integrated_data.rds",method="VIPCCA",
																						batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")

#scxx_cvae2_res <- calculate_alignment_metric(paste0(result_path,"scxx_cvae2/scmb_cvae2_3_128_ref/"),"integrated_data.rds",method="scxx_cvae2",
#																						 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
#scxx_cvae3_1_res <- calculate_alignment_metric(paste0(result_path,"scxx_cvae3/"),"output.h5ad",method="scxx_cvae3_fz",
#																							 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:32,integrated_file="integrated_processed_fz.rds")
#scxx_cvae3_2_res <- calculate_alignment_metric(paste0(result_path,"scxx_cvae3/"),"output.h5ad",method="scxx_cvae3_sz",
#																							 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx2",dims=1:32,integrated_file="integrated_processed_sz.rds")
scanorama_res <- calculate_alignment_metric(paste0(result_path, "scanorama/"),"integrated_data.rds",method="Scanorama",
																						batch_col=batch_col,celltype_col=celltype_col,reduction="scanorama",dims=1:16,fmt="rds")
desc_res <- calculate_alignment_metric(paste0(result_path, "desc/"),"integrated_data.rds",method="DESC",
																			 batch_col=batch_col,celltype_col=celltype_col,reduction="Embeded_z0.3",dims=1:16,fmt="rds")
seurat_res <- calculate_alignment_metric(paste0(result_path, "seurat/"),"seurat_v3_integrated_data.rds",method="Seurat V3",
																				 batch_col=batch_col,celltype_col=celltype_col,reduction="pca",dims=1:16,fmt="rds")
scalign_res <- calculate_alignment_metric(paste0(result_path, "scalign/"),"scalign_pancreas_integrated.rds",method="ScAlign",
																					batch_col=batch_col,celltype_col=celltype_col,reduction="ALIGNED.MultiCCA",dims=1:16,fmt="rds")
mnn_res <- calculate_alignment_metric(paste0(result_path, "mnn/"),"mnn_integrated_data.rds",method="MNN",
																			batch_col=batch_col,celltype_col=celltype_col,reduction="pca",dims=1:16,fmt="rds")
liger_res <- calculate_alignment_metric(paste0(result_path, "liger/"),"liger_integrated_data.rds",method="LIGER",
																				batch_col=batch_col,celltype_col=celltype_col,reduction="iNMF",dims=1:16,fmt="rds")
harmony_res <- calculate_alignment_metric(paste0(result_path, "harmony/"),"integrated_data.rds",method="Harmony",
																					batch_col=batch_col,celltype_col=celltype_col,reduction="harmony",dims=1:16,fmt="rds")

# plot_metric(scxx_cvae_res,scxx_cvae2_res,scxx_cvae3_1_res,scxx_cvae3_2_res,
# 						scanorama_res,desc_res,seurat_res,mnn_res, result_path=result_path)

datasets.integrated<-readRDS(paste0(result_path, "/scxx_cvae/cvae_var2_10_128_16/integrated_processed.rds"))
datasets.integrated$dataset<-
	as.factor(datasets.integrated$dataset)
levels(datasets.integrated$dataset)<-c("Celseq","Celseq2","Fluidigmc1","Indrop1","Indrop2","Indrop3","Indrop4","Smartseq2")
datasets.integrated$celltype<-
	as.factor(datasets.integrated$celltype)
levels(datasets.integrated$celltype)<-c("Acinar","Activated stellate","Alpha", "Beta", "Delta", "Ductal","Endothelial","Epsilon","Gamma", "Macrophage","Mast","Quiescent stellate","Schwann")


celltypes <- as.factor(datasets.integrated@meta.data[[celltype_col]])
batches <- datasets.integrated@meta.data[[batch_col]]
tc<-data.frame(batches=batches, celltypes=celltypes)

batch_type = unique(batches)
lz=60
gbatch<-ggplot(tc,aes(batches))+geom_bar(aes(fill=celltypes))+
	theme(legend.position = "bottom",legend.direction = "horizontal",legend.title = element_blank(),panel.grid.major = element_blank(), axis.text=element_text(size=lz,face="bold"),
				axis.text.x =element_text(size=lz,face="bold",angle = 45,vjust = 0.5),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
				legend.key.size = unit(2,"cm"),legend.text = element_text(size=lz,face = "bold"))
gcelltype<-ggplot(tc,aes(celltypes))+geom_bar(aes(fill=batches))+
	theme(legend.position = "bottom",legend.direction = "horizontal",legend.title = element_blank(),panel.grid.major = element_blank(), axis.text=element_text(size=lz,face="bold"),
				axis.text.x =element_text(size=lz,face="bold",angle = 45,vjust = 0.5),axis.title = element_blank(),legend.text = element_text(size=lz,face = "bold"),
				legend.key.size = unit(2,"cm"),panel.grid.minor = element_blank(),panel.background = element_blank())

png(paste0(result_path,"batch_composition.png"),width = 40, height = 50, units = "in", res = 100)
plot_grid(gbatch,gcelltype,ncol = 1,labels = "auto",label_size = 100,label_x = 0, label_y = 1,hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,1),"cm"))
dev.off()

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(8)
p1<-DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col) +
	scale_color_manual(values=mycolors)+
	theme(legend.text=element_text(size=70,face = "bold"),legend.position = c(0.1,0.5))+
	guides(colour = guide_legend(override.aes = list(size=15),nrow = 1))
legend1 <- cowplot::get_legend(p1)
p2<-DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col) +
	theme(legend.text=element_text(size=70,face="bold"),legend.position = c(0.05,0.5))+
	guides(colour = guide_legend(override.aes = list(size=15),nrow = 2))
legend2 <- cowplot::get_legend(p2)

ari_data<-c()
sil<-c()
mm<-c()
tt<-c()
rr<-c()
lrr<-c()
i=0
gb<-list()
gc<-list()
gd<-list()
ge<-list()
df_kbet<-list()
library(rlist)
alg_names<-c("DESC","Harmony", "LIGER","MNN","ScAlign",
						 "Scanorama","Seurat V3","VIPCCA")
for(res in list(desc_res,harmony_res,liger_res,mnn_res,
                scalign_res,scanorama_res,seurat_res,vipcca_res)){
	i=i+1
	mm <- c(mm,res[[1]])
	sil <- c(sil,res[[2]])
	ari_data <- rbind(ari_data,res[[3]]) 
	tt<- c(tt,length(res[[2]]))
	rr<-c(rr,res[[6]][[1]][,3])
	lrr<-c(lrr,length(res[[6]][[1]][,3]))
	gb<-list.append(gb,res[[4]])
	gc<-list.append(gc,res[[5]])
	gd<-list.append(gd,res[[6]][[2]])
	ge<-list.append(ge,res[[7]][[2]])
	df_kbet<-list.append(df_kbet,res[[7]][[1]])
}
df_kbet<-do.call("rbind",df_kbet)
df_kbet[,"methods"]=rep(alg_names,each=65)
df_kbet$celltype<-as.factor(df_kbet$celltype)
levels(df_kbet$celltype) <- c("Acinar","Activated stellate","Alpha", "Beta", "Delta", "Ductal","Endothelial","Epsilon","Gamma", "Macrophage","Mast","Quiescent stellate","Schwann")

# ga<-list(scxx_cvae_res[[4]],
# 				 # scxx_cvae2_res[[4]], scxx_cvae3_1_res[[4]],scxx_cvae3_2_res[[4]],
# 				 scanorama_res[[4]],desc_res[[4]],seurat_res[[4]],scalign_res[[4]],mnn_res[[4]],liger_res[[4]],scvi_res[[4]],
# 				 scxx_cvae_res[[5]],#scxx_cvae2_res[[5]],
# 				 # scxx_cvae3_1_res[[5]],scxx_cvae3_2_res[[5]],
# 				 scanorama_res[[5]],desc_res[[5]],seurat_res[[5]],scalign_res[[5]],mnn_res[[5]],liger_res[[5]],scvi_res[[5]])
# gd<-list(scxx_cvae_res[[6]][[2]],# scxx_cvae2_res[[6]][[2]],
# 				 scanorama_res[[6]][[2]],desc_res[[6]][[2]],seurat_res[[6]][[2]],
# 				 scalign_res[[6]][[2]],mnn_res[[6]][[2]],liger_res[[6]][[2]],scvi_res[[6]][[2]])


ari_data$group <- rep(c("DESC","Harmony","LIGER","MNN","ScAlign",
									 "Scanorama","Seurat V3","VIPCCA"),
								 each=10)
ari_data$ari_batch<-1-ari_data$ari_batch

# with solution=0.8 by default in FindClusters
ss<-seq(from=3, to = 80,by = 10)
ari_data<-ari_data[ss,]
sil_data<- data.frame(csil=sil,group=rep(c("DESC","Harmony","LIGER","MNN","ScAlign",
																					 "Scanorama","Seurat V3","VIPCCA"),
																				 tt))
kbet_data<- data.frame(mkb=rr,group=rep(c("DESC","Harmony","LIGER","MNN","ScAlign",
																					"Scanorama","Seurat V3","VIPCCA"),
																					lrr),
											neighbor_size=rep(seq(5,25,by = 5),8))
# plot grid
#result_path
#gb<-list.append(gb,legend1)
# png(paste0(result_path,"dataset_umap.png"),width=30, height = 10, units="in", res=100)
# g1row<-plot_grid(plotlist=gb, labels =c("a","b","c","d","e","f","g","h"),label_size = 40,nrow = 2)
# g1row<-plot_grid(g1row,legend1,rel_widths = c(4,0.5))
# plot(g1row)
# dev.off()

png(paste0(result_path,"dataset_umap.png"),width=50, height = 20, units="in", res=100)
g1row<-plot_grid(plotlist=gb, nrow = 2)
# 	g1row<-g1row+draw_label("UMAP-1",x = 0.4,y=0,hjust = 0,size = 70,fontface = "bold")+
# 	draw_label("UMAP-2",	x = 0,y=0.4,hjust = 0,size = 70,angle = 90,fontface = "bold")+
# 	theme(plot.margin = margin(0, 0, 100, 100,unit = "pt"))
g1row<-plot_grid(g1row,legend1,rel_heights =  c(6,1),ncol = 1)
plot(g1row)
dev.off()

png(paste0(result_path,"celltype_umap.png"),width=50, height = 20, units="in", res=100)
g2row<-plot_grid(plotlist=gc, nrow = 2)
# g2row<-g2row+draw_label("UMAP-1",x = 0.4,y=0,hjust = 0,size = 70,fontface = "bold")+
# 	draw_label("UMAP-2",	x = 0,y=0.4,hjust = 0,size = 70,angle = 90,fontface = "bold")+
# 	theme(plot.margin = margin(0, 0, 100, 100,unit = "pt"))
g2row<-plot_grid(g2row,legend2,rel_heights = c(6,1.0),ncol=1)
plot(g2row)
dev.off()
#png(paste0(result_path,"cell_alignment_umap.png"),width=32, height = 8, units="in", res=100)
#cp<-plot_grid(plotlist = ga, labels = "auto", nrow=2)
#cp2<-plot_grid(cp,plot_grid(legend1,legend2,ncol=1), rel_widths = c(9,0.5))
#dev.off()

# plot metric

df=data.frame(mm=mm,aa=c("DESC","Harmony","LIGER","MNN","ScAlign",
												 "Scanorama","Seurat V3","VIPCCA"))
colors=c("cyan","chocolate","coral","chartreuse",
				 "burlywood","blue","orchid","deeppink1")
g4<-ggplot(data = df,aes(x=aa,y=mm))+
	geom_bar(stat="identity",fill=colors)+
	ylab("Mixing metric")+
	theme(axis.title.x = element_blank(),
				axis.text  = element_text(size = 60,face = "bold"),
				axis.text.x  = element_blank(),
				axis.title.y = element_text(size = 70,face="bold"),
				legend.title = element_blank(),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
g5<-ggplot(data = kbet_data, aes(x=neighbor_size,y=mkb))+
	geom_line(aes(group=group,color=group),size=10)+ylab("kBET")+
	scale_color_manual(values=colors)+
	xlab("% Sample size")+
	theme(legend.title = element_blank(),legend.direction = "horizontal",
				legend.key.size = unit(4,"cm"),legend.key.height = unit(2,"cm"),
				legend.text = element_text(size = 70,face="bold"),
				axis.text = element_text(size = 60,face = "bold"),
				axis.title = element_text(size = 70,face="bold"),legend.position = c(0.5,0.5),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
	guides(colour = guide_legend(override.aes = list(size=15),nrow=1))
legend3<-cowplot::get_legend(g5)
g5<-g5+NoLegend()
g3row<-plot_grid(plotlist = list(g4,g5),
								 labels=c("c","d",""),nrow=1,label_size = 100,rel_widths = c(2,2,1),
								 label_x = 0, label_y = 1,
								 hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,3),"cm"))
png(paste0(result_path,"mm_kbet.png"),width=50, height = 15, units="in", res=100)
plot(g3row)
dev.off()



g1<-ggplot(data = ari_data, aes(y=ari_celltype,x=group))+
	geom_bar(stat="identity",fill=colors)+
	ylab("ARI")+
	theme(axis.title.x = element_blank(),
				axis.text  = element_text(size = 60,face = "bold"),
				axis.text.x  = element_blank(),
				axis.title.y = element_text(size = 70,face="bold"),
				legend.title = element_blank(),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())

# # Compute or Extract Silhouette Information
# g2<-ggplot(data = sil_data, aes(x=group,y=csil))+
# 	geom_boxplot()+
# 	ylab("Silhouette")+
# 	theme(axis.title.x = element_blank())

# # kBET line plot
# kbet_data_sub<-kbet_data[kbet_data$neighbor_size==0.5,]
# # g5<-ggplot(data= kbet_data_sub,aes(x=cell_type,y = mkb))+
# # 	geom_line(aes(color=group))+
# # 	ylab("kBET acceptance rate")+
# # 	xlab("cell types")+
# # 	theme(legend.title = element_blank())
# # legend3 <- cowplot::get_legend(g5)
# # g5<-g5+theme(legend.position = "none")
# #
# kbet_tt<-as.data.frame(table(kbet_data_sub[kbet_data_sub$mkb>0.7,]$group))
# g3<-ggplot(data = kbet_tt, aes(x=Var1,y=Freq))+
# 	geom_bar(stat = "identity",fill="orange")+
# 	ylab("# cell types")+theme_bw(base_size=30)+
# 	theme(axis.text.x = element_text(size=15))+
# 	theme(axis.title.x = element_blank())
# png(paste0(result_path,"kbet_bars.png"),width=34,height=2, units="in",res=100)
# g3
# dev.off()


# kBET test metric
png(paste0(result_path,"kbet_celltype.png"),width=28, height = 20, units="in", res=100)
plot_grid(plotlist = gd, labels = "auto", ncol = 1)
dev.off()

png(paste0(result_path,"kbet_celltype2.png"),width=35, height = 35, units="in", res=100)
plot_grid(plotlist = ge, labels = "auto", ncol = 1, label_size = 40)
dev.off()
png(paste0(result_path,"kbet_celltype3.png"),width=15, height = 20, units="in", res=100)
ggplot(df_kbet,aes(x=nb,y=ymedian,group=methods))+geom_line(aes(color=methods),size=5)+
	scale_color_manual(values=colors)+
	theme(text = element_text(size=30,face = "bold"),panel.grid.major = element_blank(),legend.title = element_blank(),
				panel.grid.minor = element_blank(),panel.background = element_blank(),legend.key.size = unit(2,"cm"),
					axis.line = element_line(colour = "black"),legend.position = "bottom",
				)+labs(x="% Sample size",y="kBET")+
	facet_wrap(~celltype,ncol = 3)
dev.off()
#
library(rlist)
#rm(list = ls())
source("./code/run/evaluation_function.R")
result_path="./results/pancreas/"
batch_col="dataset"
celltype_col="celltype"
meth_names<-c("VIPCCA","Seurat V3","Scanorama","MNN")
vipcca_adata<-readRDS(paste0(result_path,"scxx_cvae/cvae_var2_10_128_16/integrated_data.rds"))
vipcca_adata <- CreateSeuratObject(
	counts = vipcca_adata@assays$integrated@scale.data,
	meta.data = vipcca_adata@meta.data)
seurat_adata <- readRDS(paste0(result_path,"seurat/seurat_v3_integrated_data.rds"))
seurat_adata <- CreateSeuratObject(counts = seurat_adata@assays$integrated@data,
																	 meta.data = seurat_adata@meta.data)
scanorma_adata <- readRDS(paste0(result_path,"scanorama/integrated_data.rds"))
scanorma_adata <- CreateSeuratObject(counts = scanorma_adata@assays$integrated@scale.data,
																		meta.data = scanorma_adata@meta.data)
mnn_adata <- readRDS(paste0(result_path,"mnn/mnn_integrated_data.rds"))
br.all<-data.frame()
i=0
for(data.integrated in c(vipcca_adata,
												 seurat_adata,
												 scanorma_adata,
												 mnn_adata
												 )){
	i<- i+1
	VariableFeatures(data.integrated)<-rownames(data.integrated)
	for(ct in
			c("acinar","alpha","beta","ductal")
			){
		br<-run_dge_analyse(data.integrated,
												batch1 = "celseq",
												batch2 = "celseq2",
												celltype1 = ct,
												celltype_col = celltype_col,batch_col = batch_col,method = meth_names[i])
		if(!is.null(br)){
			br<-tibble::rownames_to_column(br, "genes")
			br.all<-rbind(br.all,br)
		}
		else{
			print("Cells are less than 50 in at least one group.")
		}
		print(ct)
	}
}
table(br.all$method,br.all$ctype)
saveRDS(br.all,file = "./results/pancreas/result_R/sig_genes_all.rds")
br.all<-readRDS("./results/pancreas/result_R/sig_genes_all.rds")
br.all$log_p_val<- -log10(br.all$p_val+1e-300)
g6<-ggplot(br.all, aes(y=log_p_val,x=ctype,fill=method))+
	geom_boxplot()+
	scale_fill_manual(values = colors[c(4,6,7,8)])+
	labs(y="-log (p-value)")+
	theme(legend.title = element_blank(),
				legend.text = element_blank(),
				axis.text = element_text(size = 60,face = "bold"),
				axis.title.x = element_blank(),
				axis.title = element_text(size = 70,face = "bold"),legend.position = "none",
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
	theme(legend.key.size = unit(2,"cm"))
#	guides(colour = guide_legend(override.aes = list(size=30)))


aggregate(log_p_val~method,data=br.all,median)

br.selected <- br.all[br.all$p_val_adj<0.05,]
sig_genes<-as.data.frame(table(br.selected$method,br.selected$ctype))
# sig_genes$Var1<-as.character(sig_genes$Var1)
# sig_genes$Var2<-as.character(sig_genes$Var2)
# sig_genes$Freq<-as.numeric(sig_genes$Freq)

#source("./code/run/qqplot.R")
# df_qq<-data.frame()
# ct="beta"
# for(md in levels(factor(br.all$method))){
# 	qq<-createData(md,br.all[br.all$ctype==ct & br.all$method==md,"p_val"])
# 	df_qq<-rbind(df_qq,qq)
# }
# p <- ggplot(br.all, aes(x=expected,y=observed,group=Method))
# p<-p+geom_ribbon(aes(ymin=clower,ymax=cupper),fill = "grey70")+
# 	geom_line(aes(colour=Method))
# +theme(legend.position = c(0.1,0.8))+
# 	theme(legend.title = element_blank())


# g6<-ggplot(data=sig_genes,aes(x=Var2,y=Freq))+
# 	geom_point(aes(color=Var1),size=8)+theme_bw(base_size=30)+
# 	labs(y="# significant genes",title = "celseq vs celseq2")+
# 		theme(axis.title.x = element_blank())+
# 		theme(axis.text.x = element_text(size=20))+
# 		theme(legend.title =element_blank())+
# 		theme(legend.position = c(0.3,0.6))
	

legend4<-cowplot::get_legend(g6)
g6<-g6+NoLegend()
g4row<-plot_grid(plotlist = list(g1,g6),
								 labels=c("e","f"),nrow=1,label_size = 100,rel_widths = c(2,2),
								 label_x = 0, label_y = 1,
								 hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,3),"cm"))
png(paste0(result_path,"ari_pvalue.png"),width=50, height = 15, units="in", res=100)
plot(g4row)
dev.off()


png(paste0(result_path,"all_in_one.png"),width = 50,height = 70,units = "in",res = 100)
plot_grid(plotlist = list(g1row,g2row,g3row,g4row,legend3),ncol = 1,rel_heights = c(2,2,1.5,1.5,0.2),labels = c("a","b",""),label_size = 100)+theme(plot.margin = unit(c(2,0,2,0),"cm"))
dev.off()


source("./code/run/evaluation_function.R")
br.two<-data.frame()
br2.two<-data.frame()
i<-0

for(data.integrated in c(vipcca_adata,
												 seurat_adata,
												 scanorma_adata,
												 mnn_adata
)){
	i<-i+1
	br<-run_dge_analyse(data.integrated,
											batch1 = "celseq",
											celltype1 = "alpha",
											celltype2 = "beta",
											celltype_col = celltype_col,
											batch_col = batch_col,
											method = meth_names[i])
	br2<-run_dge_analyse(data.integrated,
											batch1 = "celseq",
											batch2 = "celseq2",
											celltype1 = "alpha",
											celltype2 = "beta",
											celltype_col = celltype_col,
											batch_col = batch_col,
											method = meth_names[i])
	
	if(!is.null(br)){
		br<-tibble::rownames_to_column(br, "genes")
		br.two<-rbind(br.two,br)
	}
	else{
		print("Cells are less than 50 in at least one group.")
	}
	if(!is.null(br2)){
		br2<-tibble::rownames_to_column(br2, "genes")
		br2.two<-rbind(br2.two,br2)
	}
	else{
		print("Cells are less than 50 in at least one group.")
	}
}
saveRDS(br.two,"./results/pancreas/result_R/br_two_2.rds")
saveRDS(br2.two,"./results/pancreas/result_R/br2_two_2.rds")
br.two<-readRDS("./results/pancreas/result_R/br_two_2.rds")
br2.two<-readRDS("./results/pancreas/result_R/br2_two_2.rds")
br.two$scenarios="test2"
br2.two$scenarios="test3"
br<-rbind(br.two,br2.two)
br.selected<- br[br$p_val_adj<5e-2,]

j_index<-c()
mds<-c()
for(md in levels(factor(br.selected$method))){
	print(md)
	geneset1<-br.selected[br.selected$method==md & br.selected$scenarios=="test2","genes"]
	geneset2<-br.selected[br.selected$method==md & br.selected$scenarios=="test3","genes"]
	iset <- intersect(geneset1,geneset2)
	uset <- union(geneset1,geneset2)
	j_index<-c(j_index,(1.0*length(iset))/length(uset))
	mds<-c(mds,md)
}

br.two.selected <- br.two[br.two$p_val_adj<5e-2,]
br2.two.selected <- br2.two[br2.two$p_val_adj<5e-2,]



sig_genes_br_two<-as.data.frame(table(br.two.selected$method))
sig_genes_br2_two<-as.data.frame(table(br2.two.selected$method))


sig_genes_br_two$group="Scenario 2"
sig_genes_br2_two$group="Scenario 3"
one_df=rbind(sig_genes_br_two,sig_genes_br2_two)
g7<-ggplot(data = one_df,aes(x = Var1, y=Freq,fill=group))+
	geom_bar(position="dodge", stat="identity")+
	#geom_text(aes(label=j_index), vjust=1.6, color="white", size=3.5)+
	ylab("# significant genes")+
	theme_bw(base_size = 30)+
	theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20))+
	theme(legend.position=c(0.15,0.8),legend.title = element_blank(),
				legend.text=element_text(size=20))



g8<-ggplot(data = data.frame(sim=j_index,mds=mds),aes(x = mds, y=sim))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Jaccard Index (significant genes)")+
	theme_bw(base_size = 30)+
	theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20))
###

br.two<-readRDS("./results/pancreas/result_R/br_two_2.rds")
br2.two<-readRDS("./results/pancreas/result_R/br2_two_2.rds")
br.two <- br.two[order(br.two$p_val_adj),]
br.two.selected <- by(br.two, br.two$method, head, n=100)
br2.two <- br2.two[order(br2.two$p_val_adj),]
br2.two.selected <- by(br2.two, br2.two$method, head, n=100)
j_index <- c()
for(i in 1:4){
	geneset1<-br.two.selected[[i]][[1]]
	geneset2<-br2.two.selected[[i]][[1]]
	iset <- intersect(geneset1,geneset2)
	uset <- union(geneset1,geneset2)
	j_index<-c(j_index,(1.0*length(iset))/length(uset))
}

mds<-c("MNN","Scanorama","Seurat V3", "VIPCCA")

g9<-ggplot(data = data.frame(sim=j_index,mds=mds),aes(x = mds, y=sim))+
	geom_bar(stat="identity",fill=colors[c(4,6,7,8)])+
	ylab("Jaccard Index (significant genes)")+
	theme_bw(base_size = 30)+
	theme(axis.title.x = element_blank(),axis.text = element_text(size=40,face = "bold"),
				axis.title.y = element_text(size=40,face="bold"),
				panel.grid.major = element_blank(),panel.grid.minor = element_blank())





png(paste0(result_path,"sp_1.png"),width = 20,height = 10,units = "in",res = 100)
g9
#plot_grid(plotlist = list(g7,g8),nrow=1,labels = c("a","b"),label_size = 40)
dev.off()

###########
rm(list = ls())
library(clusterProfiler)
br.two<-readRDS("./results/pancreas/result_R/br_two_2.rds")
br2.two<-readRDS("./results/pancreas/result_R/br2_two_2.rds")
br.two.selected <- br.two[br.two$p_val_adj<5e-2,]
br2.two.selected <- br2.two[br2.two$p_val_adj<5e-2,]
library("org.Hs.eg.db")
df_kk<-data.frame()
df_kk_sp1<-data.frame()
for(mt in c("VIPCCA","MNN","Scanorama","Seurat V3")){
	geneset2<-br.two.selected[br.two.selected$method==mt,]$genes
	geneset2<-unlist(mget(x=geneset2,envir=org.Hs.egALIAS2EG,ifnotfound = NA))
	geneset2<-geneset2[!is.na(geneset2)]
	geneset3<-br2.two.selected[br2.two.selected$method==mt,]$genes
	geneset3<-unlist(mget(x=geneset3,envir=org.Hs.egALIAS2EG,ifnotfound = NA))
	geneset3<-geneset3[!is.na(geneset3)]
	
	geneset <- union(geneset2,geneset3)
	
	# kk2 <- as.data.frame(enrichKEGG(geneset2, organism="hsa", pAdjustMethod="BH", 
	# 																qvalueCutoff=0.05))
	# kk3 <- as.data.frame(enrichKEGG(geneset3, organism="hsa", pAdjustMethod="BH", 
	# 																qvalueCutoff=0.05))
	kk <- as.data.frame(enrichKEGG(geneset, organism="hsa", pAdjustMethod="BH", 
	                                qvalueCutoff=0.05))
	kk$report <- mt
	
	# kk2$report <- paste0(mt, "-s2")
	# kk3$report <- paste0(mt, "-s3")
	# df_kk<-rbind(df_kk,rbind(kk2,kk3))
	# 
	# iset <- intersect(kk2$ID,kk3$ID)
	# uset <- union(kk2$ID,kk3$ID)
	# df_kk_sp1<-rbind(df_kk_sp1,list(length(kk2$ID),length(kk3$ID),(1.0*length(iset))/length(uset), length(iset)))

	df_kk<- rbind(df_kk,kk)
}

df_kk$tqvalue <- -log10(df_kk$qvalue)
st_kk<-as.data.frame(table(df_kk$Description))
df_kk$highlight <- "conserved"
# df_kk[df_kk$Description=="Cytokine-cytokine receptor interaction","highlight"]="positive"
# df_kk[df_kk$Description=="Malaria","highlight"]="positive"
# df_kk[df_kk$Description=="Osteoclast differentiation","highlight"]="positive"
# df_kk[df_kk$Description=="Insulin secretion","highlight"]="positive"
# df_kk[df_kk$Description=="Gastric acid secretion","highlight"]="positive"
# df_kk[df_kk$Description=="cAMP signaling pathway","highlight"]="positive"

df_kk[df_kk$Description %in% as.character(st_kk[st_kk$Freq<2,"Var1"]),"highlight"]="interest"

g9 <- ggplot(df_kk, aes(y=reorder(Description,tqvalue),x=tqvalue, fill=highlight))+
	geom_bar(stat = "identity")+
	scale_x_continuous(trans = "log10")+labs(x="-log(qvalue)")+
	theme_bw(base_size = 60)+
	theme(axis.title.y = element_blank(),legend.position = "none")+
	facet_wrap(facets = .~report,nrow = 1)+
  ggtitle("Enriched KEGG pathway between alpha and beta cells in pancreatic islets")
result_path="./results/pancreas/"
png(paste0(result_path,"sp_2.png"),width = 80,height = 60,units = "in",res = 100)
g9
dev.off()


colnames(df_kk_sp1)<-c("Scenario2","Scenario3","Jaccard","overlapped")
df_kk_sp1$method <-  c("VIPCCA","MNN","Scanorama","Seurat V3")

df_kk_sp1


g10 <-ggplot(df_kk_sp1, aes(x=method,y=Jaccard))+
	geom_bar(stat = "identity",fill="orange") +
  ylab("Jaccard Index (enriched pathways)")+
	theme(axis.title = element_text(size = 30),axis.text = element_text(size = 30),
				legend.position = c(0.2,0.8),
				axis.title.x = element_blank(),
				legend.text  = element_text(size = 20),legend.title = element_text(size = 30))

png(paste0(result_path,"sp_1.png"),width = 30,height = 10,units = "in",res = 100)
plot_grid(plotlist = list(g7,g8,g10),nrow=1,labels = c("a","b","c"),label_size = 40)
dev.off()


