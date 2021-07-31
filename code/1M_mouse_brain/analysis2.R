# cpu and mem usage of vipcca and scanorama
###
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(mclust)
library(doParallel)
library(kBET)
library(FNN)
library(gridExtra)
library(grid)
source("./code/run/evaluation_function.R")
cl <- makeCluster(20)
registerDoParallel(cl)
clusterExport(cl, varlist=c("kbet_for_one_celltype", "get.knn", "kBET"), env=environment())
ann_nuclei<-readRDS("./results/1M_mouse_brain/ann_nuclei.rds")
ann_dropviz<-readRDS("./results/1M_mouse_brain/ann_dropviz.rds")
ann_sauder<-readRDS("./data/scanorama_data/data/mouse_brain/ann_saunders.rds")
ann_sauder$cluster<-unlist(lapply(strsplit(ann_sauder$subcluster,split = "-"),function(x) x[1]))
ann_sauder_simple<-ann_sauder[,c(2,3,9)]
ann_sauder_simple<-ann_sauder_simple[!duplicated(ann_sauder_simple),]
#metadata<-merge(ann_dropviz,ann_sauder,by.x="cluster",by.y="cluster",all.x=TRUE)

analyse_on_mouse_brain<-function(adata=NULL, method="VIPCCA",reduction="scxx",
																 result_path="./"){
	#browser()
	result_path<-paste0("./results/1M_mouse_brain/",method,"/", sep = "")
	if(!file.exists(paste0(result_path,method,"_um_ts_3.rds",sep = ""))){
		adata$.batch<-as.factor(adata$.batch)
		levels(adata$.batch)<-c("Nuclei","CerebellumALT",
														"Cortex noRep5 FRONTALonly","Cortex noRep5 POSTERIORonly",
														"EntoPeduncular","GlobusPallidus",
														"Hippocampus","Striatum",
														"SubstantiaNigra","Thalamus")
		adata@meta.data["tech"]="Drop-seq"
		adata$tech[adata$.batch=="Nuclei"]="SPLiT-seq"
		adata <- FindNeighbors(adata, reduction=reduction,dims = 1:16)
		adata<-RunTSNE(adata,reduction = reduction,dims=1:16,check_duplicates = FALSE)
		adata<-RunUMAP(adata,reduction=reduction,dims=1:16)
		saveRDS(adata,paste0(result_path,method,"_um_ts_3.rds",sep=""))
	}else{
		adata<-readRDS(paste0(result_path,method,"_um_ts_3.rds",sep=""))
	}
	
	# browser()
	Idents(adata)<-".batch"
	adata_splitseq<-subset(adata, idents = "Nuclei")
	adata_splitseq<-AddMetaData(adata_splitseq,metadata = ann_nuclei)
	adata_dropviz<-adata[,Idents(adata)!="Nuclei"]
	adata_dropviz<-AddMetaData(adata_dropviz, metadata = ann_dropviz)
	# cellnames<-c(colnames(adata_dropviz),colnames(adata_splitseq))
	# adata<-adata[,cellnames]
	Idents(adata)<-"tech"
	preduction="umap"
	g1<-DimPlot(adata,reduction = preduction,group.by = "tech",order = "SPLiT-seq",cols=c("grey","red"))+ggtitle(method)+
		theme_bw(base_size=30)+NoLegend()+
		theme(axis.title = element_blank(),legend.text = element_text(size = 40,face = "bold"),title = element_text(size = 70,face = "bold"),
					axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5),
					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
	g2<-DimPlot(adata,reduction = preduction,group.by = ".batch",order = "Nuclei")+ggtitle(method)+
		theme_bw(base_size=30)+NoLegend()+
		theme(axis.title = element_blank(),legend.text = element_text(size = 40,face="bold"),title = element_text(size = 70,face = "bold"),
					axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5),
					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
	g3<-DimPlot(adata_dropviz,reduction = preduction,group.by = "cluster")+
		ggtitle(paste(method,"drop-seq"))+theme_bw(base_size=30)+NoLegend()+
		theme(axis.title = element_blank(),legend.text = element_text(size = 40),title = element_text(size = 70),
					axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5),
					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
	g4<-DimPlot(adata_splitseq,reduction = preduction,group.by = "cell_type")+
		ggtitle(paste(method,"SPLiT-seq"))+theme_bw(base_size=30)+NoLegend()+
		theme(axis.title = element_blank(),legend.text = element_text(size = 40),title = element_text(size = 70),
					axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5),
					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
	
	adata_dropviz_NA_removed<-adata_dropviz[,!is.na(adata_dropviz$cluster)]
	barcode_dropviz <- colnames(adata_dropviz_NA_removed)
	barcode_splitseq<- colnames(adata_splitseq)
	celltype_dropviz <- ann_dropviz[barcode_dropviz,1]
	celltype_splitseq<- ann_nuclei[barcode_splitseq,2]
	
	Idents(adata)<-".batch"
	if(!file.exists(paste(result_path,method,"_ari_2.rds",sep = ""))){
		ari_score<-c()
		for(resolution in seq(0.1,1,0.1)){
			adata<-FindClusters(adata,resolution = resolution)
			cluster_dropviz <- adata@meta.data[barcode_dropviz,"seurat_clusters"]
			cluster_splitseq <- adata@meta.data[barcode_splitseq,"seurat_clusters"]
			ari_score_dropviz<-adjustedRandIndex(celltype_dropviz,cluster_dropviz)
			ari_score_splitseq<-adjustedRandIndex(celltype_splitseq,cluster_splitseq)
			ari_score<-c(ari_score, ari_score_dropviz,ari_score_splitseq)
		}
		df_ari<-data.frame(matrix(ari_score,ncol=2,byrow=TRUE))
		colnames(df_ari)<-c("dropseq","splitseq")
		df_ari$method <- method
		saveRDS(df_ari,paste0(result_path,method,"_ari_2.rds",sep=""))
	}else{
		df_ari<-readRDS(paste0(result_path,method,"_ari_2.rds",sep=""))
	}
	
	adata@meta.data[["is_selected"]]="all_cells_types"
	kfile=paste0("kbet_",method,"_4.rds",sep="")
	celltype_col="is_selected"
	batch_col="tech"
	if(!file.exists(paste0(result_path,kfile,sep=""))){
		require('FNN')
		com=run_kbet_for_each_celltype4(adata,
																		reduction = reduction,
																		celltype_col = celltype_col,
																		batch_col = batch_col,
																		method = method,
																		downsample = TRUE)
		saveRDS(com, paste0(result_path,kfile,sep=""))
	}else
	{
		com <- readRDS(paste0(result_path,kfile,sep=""))
	}
	g5<-com[[2]]+ylim(c(0,0.5))+theme_bw(base_size=30)+theme(axis.text.x = element_text(size=15))+
		theme(axis.title.x = element_blank())+ylab("kBET")
	
	#browser()
	max.k=300
	if(!file.exists(paste(result_path,method,"_mm_2.rds",sep = ""))){
		mm <- max.k-MixingMetric(adata, grouping.var=batch_col, 
													 reduction=reduction, dims=1:16, max.k=max.k)
		saveRDS(mm,paste(result_path,method,"_mm_2.rds",sep = ""))
	}else{
		mm<-readRDS(paste(result_path,method,"_mm_2.rds",sep = ""))
	}
	
	return(list(g1,g2,g3,g4,com,median(mm),df_ari))
}

analyse_on_mouse_brain_legend<-function(adata=NULL, method="VIPCCA",reduction="scxx",
																 result_path="./"){
	# browser()
	result_path<-paste0("./results/1M_mouse_brain/",method,"/", sep = "")
	adata<-readRDS(paste0(result_path,method,"_um_ts_3.rds",sep=""))
	#browser()
	Idents(adata)<-".batch"
	adata_splitseq<-subset(adata, idents = "Nuclei")
	adata_splitseq<-AddMetaData(adata_splitseq,metadata = ann_nuclei)
	adata_dropviz<-adata[,Idents(adata)!="Nuclei"]
	adata_dropviz<-AddMetaData(adata_dropviz, metadata = ann_dropviz)

	Idents(adata)<-"tech"
	preduction="umap"
	g1<-DimPlot(adata,reduction = preduction,group.by = "tech",order = "SPLiT-seq",cols=c("grey","red"))+ggtitle(method)+
		theme_bw(base_size=30)+
		theme(legend.text = element_text(size=70,face="bold"),legend.position = c(0.5,0.5),legend.direction = "horizontal")+
		guides(colour = guide_legend(override.aes = list(size=15),nrow = 1))
	g1<-cowplot::get_legend(g1)
	g2<-DimPlot(adata,reduction = preduction,group.by = ".batch")+ggtitle(method)+
		theme_bw(base_size=30)+theme(legend.text = element_text(size=70,face = "bold"),
																 legend.position = c(0.5,0.5),legend.direction = "horizontal")+
		guides(colour = guide_legend(override.aes = list(size=15),nrow = 2))
	g2<-cowplot::get_legend(g2)
	g3<-DimPlot(adata_dropviz,reduction = preduction,group.by = "cluster")+
		ggtitle(paste(method,"drop-seq"))+
		theme(legend.text = element_text(size=70,face = "bold"),legend.direction = "horizontal",
					legend.position = c(0.2,0.5))+
		guides(colour = guide_legend(override.aes = list(size=15),nrow = 1))
	g3<-cowplot::get_legend(g3)
	g4<-DimPlot(adata_splitseq,reduction = preduction,group.by = "cell_type")+
		ggtitle(paste(method,"SPLiT-seq"))+
		theme(legend.text = element_text(size=70,face = "bold"),legend.direction = "horizontal",
					legend.position = c(0.1,0.5))+
		guides(colour = guide_legend(override.aes = list(size=15),nrow = 20))
	g4<-cowplot::get_legend(g4)
	lz=20
	#reduction="umap"
	g5<-FeaturePlot(adata,reduction = preduction,features = c("MBP","ALDH1L1","PDGFRA","DOCK2"),split.by = "tech",by.col = F,combine = T) &
		theme(axis.title.x = element_blank(),axis.title.y.left = element_blank(),panel.border = element_blank(),
					axis.text = element_blank(),axis.ticks = element_blank(),title = element_text(size = lz),axis.title.y.right = element_text(size = lz))
	g6<-FeaturePlot(adata,reduction = preduction,features = c("RGS5","COL1A2","DNAH11","GLI2"),split.by = "tech",by.col = F,combine = T) &
		theme(axis.title.x = element_blank(),axis.title.y.left = element_blank(),panel.border = element_blank(),
					axis.text = element_blank(),axis.ticks = element_blank(),
					title = element_text(size = lz),axis.title.y.right = element_text(size = lz))
	g7<-FeaturePlot(adata,reduction = preduction,features = c("DEPTOR","MYBPC1","MEG3","FN1"),split.by = "tech",by.col = F,combine=T) &
		theme(axis.title.x = element_blank(),axis.title.y.left = element_blank(),panel.border = element_blank(),
					axis.text = element_blank(),axis.ticks = element_blank(),
					title = element_text(size = lz),axis.title.y.right = element_text(size = lz))
	g8<-FeaturePlot(adata,reduction = preduction,features = c("GAD1","GAD2","C1QB","TNR"),split.by="tech",by.col = F,combine = T)&
		theme(axis.title.x = element_blank(),axis.title.y.left = element_blank(),panel.border = element_blank(),
					axis.text = element_blank(),axis.ticks = element_blank(),
					title = element_text(size = lz),axis.title.y.right = element_text(size = lz))
	g9<-FeaturePlot(adata,reduction = preduction,features = c("CCDC153","GJA1","DCN","TRF"),split.by="tech",by.col = F,combine = T)&
		theme(axis.title.x = element_blank(),axis.title.y.left = element_blank(),panel.border = element_blank(),
					axis.text = element_blank(),axis.ticks = element_blank(),
					title = element_text(size = lz),axis.title.y.right = element_text(size = lz))
	g10<-FeaturePlot(adata,reduction = preduction,features = c("FLT1","SOX4","MKI67","RORB"),split.by = "tech",by.col = F,combine = T)&
		theme(axis.title.x = element_blank(),axis.title.y.left = element_blank(),panel.border = element_blank(),
					axis.text = element_blank(),axis.ticks = element_blank(),
					title = element_text(size = lz),axis.title.y.right = element_text(size = lz))
	
	# browser()
	# gene_names<-c("Mbp","Aldh1l1","Pdgfra","Dock2")
	# g5[1:4]<-lapply(1:length(gene_names), function(x) { g5[[x]] + labs(title=gene_names[x]) })
	# # g5 <- CombinePlots( g5 )
	# gene_names<-c("Rgs5","Col1a2","Dnah11","Gli2")
	# g6<-lapply(1:length(gene_names), function(x) { g6[[x]] + labs(title=gene_names[x]) })
	# # g6 <- CombinePlots( g6 )
	# gene_names<-c("Deptor","Mybpc1","Meg3","Fn1")
	# g7<-lapply(1:length(gene_names), function(x) { g7[[x]] + labs(title=gene_names[x]) })
	# # g7 <- CombinePlots( g7 )
	# gene_names<-c("Gad1","Gad2","C1qb","Tnr")
	# g8<-lapply(1:length(gene_names), function(x) { g8[[x]] + labs(title=gene_names[x]) })
	# # g8 <- CombinePlots( g8 )
	# gene_names<-c("Ccdc153","Gja1","Dcn","Trf")
	# g9<-lapply(1:length(gene_names), function(x) { g9[[x]] + labs(title=gene_names[x]) })
	# # g9 <- CombinePlots( g9 )
	# gene_names<-c("Flt1","Sox4","Mki67","Rorb")
	# g10<-lapply(1:length(gene_names), function(x) { g10[[x]] + labs(title=gene_names[x]) })
	# # g10 <- CombinePlots( g10 )

	# 								
	# g5<-FeaturePlot(adata,reduction = "tsne", features = c("MBP","ALDH1L1","PDGFRA","DOCK2",
	# 																		"RGS5","COL1A2","DNAH11","GLI2",
	# 																		"DEPTOR","MYBPC1","MEG3","FN1",
	# 																		"GAD1","GAD2","C1QB","TNR",
	# 																		"CCDC153","GJA1","DCN","TRF",
	# 																		"FLT1","SOX4","MKI67"),by.col = T,ncol = 4)

	g11<-(g5/g6)/g7
	g12<-(g8/g9)/g10
	return(list(g1,g2,g3,g4,g11,g12))
}

adata_vipcca<-readRDS("./results/1M_mouse_brain/VIPCCA/output_unq_subset_2.rds")
adata_scanorama<-readRDS("./results/1M_mouse_brain/Scanorama/output_unq_subset_2.rds")
# adata_scvi <- ReadH5AD("./results/1M_mouse_brain/scVI/output_sub_2.h5ad")
adata_desc <-  readRDS("./results/1M_mouse_brain/DESC/output_unq_subset_2.rds")
adata_harmony <- readRDS("./results/1M_mouse_brain/Harmony/output_unq_subset_2.rds")

glist_vipcca<-analyse_on_mouse_brain(adata = adata_vipcca, reduction = "scxx", method="VIPCCA")
glist_scanorama<-analyse_on_mouse_brain(adata = adata_scanorama, method = "Scanorama", reduction = "scanorama")
#glist_scvi <- analyse_on_mouse_brain(adata = adata_scvi,	method="scVI", reduction = "scvi")
glist_desc <- analyse_on_mouse_brain(adata = adata_desc,	method="DESC", reduction = "Embeded_z0.3")
glist_harmony <- analyse_on_mouse_brain(adata = adata_harmony,	method="Harmony",
																				reduction = "harmony")
glist_legends <- analyse_on_mouse_brain_legend()


library(rlist)
colors=c("cyan","chocolate","coral","chartreuse",
				 "burlywood","blue","orchid","deeppink1")

df_kbet<-rbind(glist_vipcca[[5]][[1]],glist_scanorama[[5]][[1]],glist_desc[[5]][[1]],glist_harmony[[5]][[1]])
df_kbet$method=rep(c("VIPCCA","Scanorama","DESC","Harmony"),each=5)

p0 <- ggplot(df_kbet,aes(x=nb,y=ymedian))+
	geom_line(aes(group=method,color=method),size=10)+ylab("kBET")+
	xlab("% Sample size")+ylim(0,0.5)+
	scale_color_manual(values=colors[c(1,2,6,8)])+
	theme(legend.title = element_blank(),legend.text = element_text(size = 70,face="bold"),
				panel.grid.major=element_blank(),axis.text = element_text(size=70,face="bold"),axis.title = element_text(size=70,face = "bold"),
				panel.grid.minor=element_blank(),panel.background = element_blank(),legend.key.size = unit(3,"cm"),legend.direction = "horizontal")+
	guides(colour = guide_legend(override.aes = list(size=15),nrow = 4))
legend_mix<-cowplot::get_legend(p0)
p0 <- p0+NoLegend()

df<-data.frame(mm=c(glist_desc[[6]],glist_harmony[[6]],glist_scanorama[[6]],glist_vipcca[[6]]),
							 method=c("DESC","Harmony","Scanorama","VIPCCA"))
p1<-ggplot(data = df,aes(x=method,y=mm))+
	geom_bar(stat="identity",fill=colors[c(1,2,6,8)])+
	ylab("Mixing metric")+theme_bw(base_size=70)+
	theme(axis.text.x= element_text(size=70,face = "bold",angle = 45,vjust = 0.5),axis.text.y=element_text(size=70,face="bold"),axis.ticks.x=element_blank(),axis.title = element_text(size=70,face = "bold"))+
	theme(axis.title.x = element_blank(),panel.grid.major=element_blank(),axis.text.x = element_blank(),
				panel.grid.minor=element_blank(),panel.border = element_blank())

df_ari<-do.call(rbind, list(glist_desc[[7]],glist_harmony[[7]],glist_scanorama[[7]],glist_vipcca[[7]]))
df_ari$sum <- (df_ari$dropseq+df_ari$splitseq)/2
df_ari<-df_ari[seq(2,40,10),]
p2<-ggplot(data = df_ari,aes(x=method,y=sum))+
	geom_bar(stat="identity",fill=colors[c(1,2,6,8)])+
	ylab("ARI")+theme_bw(base_size=70)+
	theme(axis.text.x = element_text(size=70,face = "bold",angle = 45,vjust = 0.5),axis.text.y=element_text(size=70,face="bold"),axis.ticks.x=element_blank(),axis.title = element_text(size=70,face = "bold"))+
	theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),axis.text.x = element_blank(),
				panel.grid.minor = element_blank(),panel.border = element_blank())

f_vip <- paste("./results/1M_mouse_brain/cpu_mem_profile_vipcca.txt", sep = "")
cnames<-c("cpu","mem")
x1<-read.delim(f_vip, header = F, sep = ",", comment.char = "#",col.names = cnames,stringsAsFactors = T)
x1$method="VIPCCA"
x1$ctime=(0:(dim(x1)[1]-1))*1.0/3600
x1$mem_In_G<-x1$mem*2.56
f_scanorama <- paste("./results/1M_mouse_brain/cpu_mem_profile_scanorama.txt", sep = "")
x2<-read.delim(f_scanorama, header = F, sep = ",", comment.char = "#",col.names = cnames,stringsAsFactors = T)
x2$method="Scanorama"
x2$ctime=(0:(dim(x2)[1]-1))*1.0/3600
x2$mem_In_G<-x2$mem*2.56

f_harmony <- paste("./results/1M_mouse_brain/cpu_mem_profile_harmony.txt", sep = "")
x3<-read.delim(f_harmony, header = F, sep = ",", comment.char = "#",col.names = cnames,stringsAsFactors = T)
x3$method="Harmony"
x3$ctime=(0:(dim(x3)[1]-1))*1.0/3600
x3$mem_In_G<-x3$mem*2.56

f_desc <- paste("./results/1M_mouse_brain/cpu_mem_profile_desc_new.txt", sep = "")
x4<-read.delim(f_desc, header = F, sep = ",", comment.char = "#",col.names = cnames,stringsAsFactors = T)
x4$method="DESC"
x4$ctime=(0:(dim(x4)[1]-1))*1.0/3600
x4$mem_In_G<-x4$mem*2.56
#
df_scale <- do.call(rbind, list(x1,x2,x3,x4))

df_cpu_mean<-aggregate(cpu~method,df_scale,FUN = mean)
df_cpu_peak<-aggregate(cpu~method,df_scale,FUN = max)
df_mem_mean<-aggregate(mem_In_G~method,df_scale,FUN = mean)
df_mem_peak<-aggregate(mem_In_G~method,df_scale,FUN = max)
df_time_used<-aggregate(ctime~method,df_scale,FUN = max)

df_cpu_mean$cputime<-df_cpu_mean$cpu*df_time_used$ctime/100

p3<-ggplot(data = df_cpu_mean,aes(x=method,y=cputime))+
	geom_bar(stat="identity",fill=colors[c(1,2,6,8)])+
	ylab("Total cputime (H)")+theme_bw(base_size=60)+
	theme(axis.text.x = element_text(size=70,face = "bold",angle = 45,vjust = 0.5))+
	theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_blank(),
				panel.border = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=70,face = "bold"),axis.text.y =element_text(size=70,face = "bold") )
p4<-ggplot(data = df_cpu_peak,aes(x=method,y=cpu))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Peak cpu (%)")+theme_bw(base_size=50)+
	theme(axis.text.x = element_text(size=30))+
	theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5<-ggplot(data = df_mem_mean,aes(x=method,y=mem_In_G))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Ave mem (G)")+theme_bw(base_size=50)+
	theme(axis.text.x = element_text(size=30))+
	theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p6<-ggplot(data = df_mem_peak,aes(x=method,y=mem_In_G))+
	geom_bar(stat="identity",fill=colors[c(1,2,6,8)])+
	ylab("Peak mem (G)")+theme_bw(base_size=70)+
	theme(axis.text.x = element_text(size=70,face = "bold",angle = 45,vjust = 0.5))+
	theme(axis.title.x = element_blank(),axis.text.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
				panel.border = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=70,face = "bold"),axis.text.y =element_text(size=70,face = "bold"))
p7<-ggplot(data = df_time_used,aes(x=method,y=ctime))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Time elapsed (H)")+theme_bw(base_size=50)+
	theme(axis.text.x = element_text(size=30))+
	theme(axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

glist_all<-c(glist_desc[c(1,2)],glist_harmony[c(1,2)],glist_scanorama[c(1,2)],glist_vipcca[c(1,2)])
glist_all<-glist_all[c(seq(1,8,2),seq(2,8,2))]
# glist_all <- lapply(glist_all, function(x) x+
# 											theme(title = element_text(size = 30),
# 														panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
g1<-plot_grid(plotlist=glist_all[1:4],ncol=4)
g1<-plot_grid(g1,glist_legends[[1]],rel_heights  = c(4,1),ncol = 1)
png("./results/1M_mouse_brain/row1.png",width = 50,height = 10,units = "in",res = 100)
plot(g1)
dev.off()
g2<-plot_grid(plotlist = glist_all[5:8],ncol=4)
g2<-plot_grid(g2,glist_legends[[2]],rel_heights  = c(4,1),ncol = 1)
png("./results/1M_mouse_brain/row2.png",width = 50,height = 10,units = "in",res = 100)
plot(g2)
dev.off()

g3<-plot_grid(p0,p1,p2,p3,nrow=1,rel_widths = c(2,2,2,2),labels=c("c","d","e","f"),label_size=100,
							label_x = 0, label_y = 1,
							hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,2),"cm"))
png("./results/1M_mouse_brain/row3.png",width = 50,height = 15,units = "in",res = 100)
plot(g3)
dev.off()

g4<-plot_grid(p6,legend_mix,NULL,nrow=1,rel_widths=c(1,1,2),label_size=100,labels=c("g"),
							label_x = 0, label_y = 1,
							hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,2),"cm"))
png("./results/1M_mouse_brain/row4.png",width = 50,height = 15,units = "in",res = 100)
plot(g4)
dev.off()

g5<-plot_grid(g1,g2,g3,g4,ncol=1,labels = c("a","b"),label_size = 100,rel_heights = c(10,10,10,10))
png("./results/1M_mouse_brain/all-in-one.png",width=50,height=60,units = "in",res=100)
g5
dev.off()

# plot_main<-plot_grid(plotlist = glist_all,ncol = 4,labels = "auto",label_size = 30)
# glist_legends_main<-list.append(glist_legends[1:2],legend_mix,NULL)
# plot_main_legend<-plot_grid(plotlist = glist_legends_main, ncol = 1)
# 
# png("./results/1M_mouse_brain/plot_all_umap.png",width = 30, height = 40,units = "in", res = 100)
# #plot_main+plot_main_legend+plot_layout(ncol = 2,widths = c(7,2))
# g4
# dev.off()

# # for tsne only
# glist_all<-c(glist_desc[c(1,2)],glist_scanorama[c(1,2)],glist_scvi[c(1,2)],glist_vipcca[c(1,2)])
# glist_all<-glist_all[c(seq(1,8,2),seq(2,8,2))]
# plot_main<-plot_grid(plotlist = glist_all,nrow = 2,labels = "auto",label_size = 30)
# glist_legends_main<-list.append(glist_legends[1:2])
# plot_main_legend<-plot_grid(plotlist = glist_legends_main, ncol = 1)
# png("./results/1M_mouse_brain/plot_all_tsne.png",width = 35, height = 10,units = "in", res = 300)
# plot_main+plot_main_legend+plot_layout(ncol = 2,widths = c(7,2))
# dev.off()

glist_all<-c(glist_desc[c(3,4)],glist_harmony[c(3,4)],glist_scanorama[c(3,4)],glist_vipcca[c(3,4)])
glist_all<-glist_all[c(seq(1,8,2),seq(2,8,2))]
glist_all <- lapply(glist_all, function(x) x+
											theme(title = element_text(size = 70,face = "bold"),
														legend.title = element_text(size = 70,face="bold"),
														panel.grid.major = element_blank(),panel.grid.minor = element_blank()))

row1<-{plot_grid(plotlist = glist_all[1:4],nrow = 1)/glist_legends[[3]]}+plot_layout(ncol = 1,heights = c(3,1))
row2<-plot_grid(plotlist = glist_all[5:8],nrow = 1)/glist_legends[[4]]+plot_layout(ncol = 1,heights = c(3,5))

# png("./results/1M_mouse_brain/sp_celltypes.png",width = 50,height = 30,units = "in",res = 100)
# plot_grid(plotlist = glist_all,nrow = 2,labels = "auto",label_size = 100)/plot_grid(plotlist = glist_legends[3:4],nrow = 1,rel_widths  = c(1:6))
# dev.off()
png("./results/1M_mouse_brain/sp_celltypes.png",width = 60,height = 50,units = "in",res = 20)
plot_grid(row1,row2,ncol = 1,rel_heights = c(4,8),labels = "auto",label_size = 100)
dev.off()

png("./results/1M_mouse_brain/sp_markers_umap_1.png",width = 10,height = 15,units = "in",res = 100)
glist_legends[[5]]
dev.off()
png("./results/1M_mouse_brain/sp_markers_umap_2.png",width = 10,height = 15,units = "in",res = 100)
glist_legends[[6]]
dev.off()


# 
# library(rlist)
# methods <- c("DESC", "Scanorama", "scVI", "VIPCCA")
# sp_list <- list()
# for(method in methods){
# 	adata<-readRDS(paste0("./results/1M_mouse_brain/",method,"/",method,"_umap_2.rds",sep=""))
# 	levels(adata@meta.data$.batch)[1]="CNS (SPLiT-seq)"
# 	adata<-RunTSNE(adata,reduction = "scanorama",dims=1:16,check_duplicates = FALSE)
# 	sp <- DimPlot(adata,reduction = "umap",group.by = ".batch")+ggtitle(method)+
# 		theme_bw(base_size=30)+theme(axis.text.x = element_text(size=15))+
# 		theme(axis.title.x = element_blank())
# 	legend_region<-cowplot::get_legend(sp)
# 	sp <- sp + NoLegend()
# 	sp_list<-list.append(sp_list,sp)
# }
# glist_all<-c(glist_desc[1:3],glist_scanorama[1:3],glist_scvi[1:3],glist_vipcca[1:3])
# glist_all<-glist_all[c(seq(1,12,3),seq(2,12,3),seq(3,12,3))]
# glist_all<-c(glist_all,sp_list)
# library(patchwork)
# png("./results/1M_mouse_brain/sp_2.png",width = 24, height = 24, units = "in", res = 300)
# plot_grid(plotlist = glist_all, ncol = 4，label="auto")
# dev.off()
# png("./results/1M_mouse_brain/sp_3.png",width = 16,height = 8,units = "in",res = 300)
# sp_list[[2]]+legend_region+plot_layout(widths=c(1,1))
# dev.off()

