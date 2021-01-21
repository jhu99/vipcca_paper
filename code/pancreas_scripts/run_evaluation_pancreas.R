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
res_list <- list()
ari<-c()
gari<-c()
rr<-c()
lrr<-c()
gd<-list()
curr_path=paste0(result_path,"scxx_cvae_1e4/")
if(!file.exists(paste0(curr_path,"kbet_cvae.png"))){
	methods<-dir(curr_path,pattern = "var2")
	for(dir_cvae in methods){
		print(dir_cvae)
		h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(curr_path,dir_cvae,"/"),
									 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
									 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
		cvae_res <- calculate_alignment_metric(paste0(curr_path,dir_cvae,"/"),"integrated_data.rds",method=dir_cvae,algorithm="cover_tree",
																								batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
		ari <- c(ari,cvae_res[[3]])
		gari <- c(gari,length(cvae_res[[3]]))
		rr<-c(rr,cvae_res[[6]][[1]][,3])
		lrr<-c(lrr,length(cvae_res[[6]][[1]][,3]))
		gd[[dir_cvae]] <- cvae_res[[6]][[2]]
	}
	ari_data_cvae<- data.frame(cari=ari,group=rep(methods, each=10))
	kbet_data_cvae<- data.frame(mkb=rr,group=rep(methods,lrr))
	png(paste0(curr_path,"kbet_cvae.png"),width=12, height = 12, units="in", res=100)
	boxplot(mkb~group, data=kbet_data_cvae, ylab="kBET (acceptance rate)", las=2)
	dev.off()
	png(paste0(curr_path,"ari_cvae.png"),width=12, height = 12, units="in", res=100)
	boxplot(cari~group,data = ari_data_cvae, ylab="Adjusted Rand Index", las=2)
	dev.off()
}

saveRDS(list(ari_data_cvae,kbet_data_cvae),file = "./results/pancreas/scxx_cvae_1e4/cvae_stat.rds")

# head(kbet_data_cvae)
# head(ari_data_cvae)
# saveRDS(list(kbet_data_cvae,ari_data_cvae),file = paste0(curr_path,"cvae_stat_1235.rds"))
#rm(list = ls())
res_list <- list()
ari<-c()
gari<-c()
rr<-c()
lrr<-c()
gd<-list()
curr_path=paste0(result_path,"scxx_cvae2/")
if(!file.exists(paste0(curr_path,"kbet_cvae2.png"))){
	methods<-dir(curr_path,pattern = "cvae2_2")
	for(dir_cvae in methods){
		print(dir_cvae)
		h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(curr_path,dir_cvae,"/"),
									 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
									 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
		cvae_res <- calculate_alignment_metric(paste0(curr_path,dir_cvae,"/"),"integrated_data.rds",method=dir_cvae,algorithm="cover_tree",
																					 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
		ari <- c(ari,cvae_res[[3]])
		gari <- c(gari,length(cvae_res[[3]]))
		rr<-c(rr,cvae_res[[6]][[1]][,3])
		lrr<-c(lrr,length(cvae_res[[6]][[1]][,3]))
		gd[[dir_cvae]] <- cvae_res[[6]][[2]]
	}
	ari_data_cvae2<- data.frame(cari=ari,group=rep(methods, each=10))
	kbet_data_cvae2<- data.frame(mkb=rr,group=rep(methods,lrr))
	png(paste0(curr_path,"kbet_cvae.png"),width=12, height = 12, units="in", res=100)
	boxplot(mkb~group, data=kbet_data_cvae2, ylab="kBET (acceptance rate)", las=2)
	dev.off()
	png(paste0(curr_path,"ari_cvae.png"),width=12, height = 12, units="in", res=100)
	boxplot(cari~group,data = ari_data_cvae2, ylab="Adjusted Rand Index", las=2)
	dev.off()
}
saveRDS(list(ari_data_cvae2,kbet_data_cvae2),file = "./results/pancreas/scxx_cvae2/cvae2_stat_1510.rds")
ax1<-readRDS("./results/pancreas/scxx_cvae/cvae_stat.rds")
ax2<-readRDS("./results/pancreas/scxx_cvae2/cvae2_stat_1510.rds")
ax3<-readRDS("./results/pancreas/scxx_cvae_1e4/cvae_stat.rds")
ad1<-aggregate(cari~group,ax1[[1]],max)
md1<-aggregate(mkb~group,ax1[[2]],median)
md1[md1$mkb>0.9,]
ad1[rownames(md1[md1$mkb>0.9,]),]

ad2<-aggregate(cari~group,ax2[[1]],max)
md2<-aggregate(mkb~group,ax2[[2]],median)
md2[md2$mkb>0.8,]
ad2[rownames(md2[md2$mkb>0.8,]),]

ad3<-aggregate(cari~group,ax3[[1]],mean)
md3<-aggregate(mkb~group,ax3[[2]],median)
md3[md3$mkb>0.9,]
ad3[rownames(md3[md3$mkb>0.9,]),]


#saveRDS(res, "./results/pancreas/scxx_cvae_res.rds")
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scxx_cvae/scxx_cvae_l1_b128/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
# h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scxx_cvae2/",
# 							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
# 							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scanorama/",
							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/desc/",
							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
h5ad_to_seurat(input_file = "output.h5ad", result_path = "./results/pancreas/scvi/",
							 annotation_file="./data/panc8/panc8_10x_mtx/barcodes.tsv",
							 col.names = c("barcode", "ident","nCount","nFeature","tech","replicate","assigned_cluster","celltype","dataset"))
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
scalign_res <- calculate_alignment_metric(paste0(result_path, "scalign/"),"scalign_pancreas_integrated.rds",method="scAlign",
																					batch_col=batch_col,celltype_col=celltype_col,reduction="ALIGNED.MultiCCA",dims=1:16,fmt="rds")
mnn_res <- calculate_alignment_metric(paste0(result_path, "mnn/"),"mnn_integrated_data.rds",method="MNN",
																			batch_col=batch_col,celltype_col=celltype_col,reduction="pca",dims=1:16,fmt="rds")
liger_res <- calculate_alignment_metric(paste0(result_path, "liger/"),"liger_integrated_data.rds",method="LIGER",
																				batch_col=batch_col,celltype_col=celltype_col,reduction="iNMF",dims=1:16,fmt="rds")
scvi_res <- calculate_alignment_metric(paste0(result_path, "scvi/"),"integrated_data.rds",method="scVI",
																				batch_col=batch_col,celltype_col=celltype_col,reduction="scvi",dims=1:16,fmt="rds")
harmony_res <- calculate_alignment_metric(paste0(result_path, "harmony/"),"integrated_data.rds",method="Harmony",
																					batch_col=batch_col,celltype_col=celltype_col,reduction="harmony",dims=1:16,fmt="rds")

# plot_metric(scxx_cvae_res,scxx_cvae2_res,scxx_cvae3_1_res,scxx_cvae3_2_res,
# 						scanorama_res,desc_res,seurat_res,mnn_res, result_path=result_path)

datasets.integrated<-readRDS(paste0(result_path, "/scxx_cvae/cvae_var2_10_128_16/integrated_processed.rds"))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(8)
p1<-DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col) +
	scale_color_manual(values=mycolors)+
	theme(legend.text=element_text(size=40))+
	guides(colour = guide_legend(override.aes = list(size=10)))
legend1 <- cowplot::get_legend(p1)
p2<-DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col) +
	theme(legend.text=element_text(size=40))+
	guides(colour = guide_legend(override.aes = list(size=10)))
legend2 <- cowplot::get_legend(p2)

ari<-c()
sil<-c()
mm<-c()
tt<-c()
rr<-c()
lrr<-c()
i=0
gb<-list()
gc<-list()
gd<-list()
library(rlist)
alg_names<-c("DESC","LIGER","MNN","scAlign",
						 "Scanorama","scVI","Seurat V3","VIPCCA")
for(res in list(desc_res,liger_res,mnn_res,scalign_res,
								scanorama_res,scvi_res,seurat_res,vipcca_res)){
	i=i+1
	mm <- c(mm,res[[1]])
	sil <- c(sil,res[[2]])
	ari <- c(ari,res[[3]])
	tt<- c(tt,length(res[[2]]))
	rr<-c(rr,res[[6]][[1]][,3])
	lrr<-c(lrr,length(res[[6]][[1]][,3]))
	gb<-list.append(gb,res[[4]])
	gc<-list.append(gc,res[[5]])
	gd<-list.append(gd,res[[6]][[2]])
}

# ga<-list(scxx_cvae_res[[4]],
# 				 # scxx_cvae2_res[[4]], scxx_cvae3_1_res[[4]],scxx_cvae3_2_res[[4]],
# 				 scanorama_res[[4]],desc_res[[4]],seurat_res[[4]],scalign_res[[4]],mnn_res[[4]],liger_res[[4]],scvi_res[[4]],
# 				 scxx_cvae_res[[5]],#scxx_cvae2_res[[5]],
# 				 # scxx_cvae3_1_res[[5]],scxx_cvae3_2_res[[5]],
# 				 scanorama_res[[5]],desc_res[[5]],seurat_res[[5]],scalign_res[[5]],mnn_res[[5]],liger_res[[5]],scvi_res[[5]])
# gd<-list(scxx_cvae_res[[6]][[2]],# scxx_cvae2_res[[6]][[2]],
# 				 scanorama_res[[6]][[2]],desc_res[[6]][[2]],seurat_res[[6]][[2]],
# 				 scalign_res[[6]][[2]],mnn_res[[6]][[2]],liger_res[[6]][[2]],scvi_res[[6]][[2]])


ari_data<- data.frame(cari=ari,group=rep(c("DESC","LIGER","MNN","scAlign",
																					 "Scanorama","scVI","Seurat V3","VIPCCA"),
																				 each=10))
# with solution=0.8 by default in FindClusters
ss<-seq(from=3, to = 80,by = 10)
ari_data<-ari_data[ss,]
sil_data<- data.frame(csil=sil,group=rep(c("DESC","LIGER","MNN","scAlign",
																					 "Scanorama","scVI","Seurat V3","VIPCCA"),
																				 tt))
kbet_data<- data.frame(mkb=rr,group=rep(c("DESC","LIGER","MNN","scAlign",
																					"Scanorama","scVI","Seurat V3","VIPCCA"),
																					lrr),
											neighbor_size=rep(seq(0.15,0.95,by = 0.05),13*8),
											cell_type=rep(rep(1:13,each=17),8))
# plot grid
#result_path
#gb<-list.append(gb,legend1)
png(paste0(result_path,"dataset_umap.png"),width=30, height = 10, units="in", res=100)
g1row<-plot_grid(plotlist=gb, labels =c("a","b","c","d","e","f","g","h"),label_size = 40,nrow = 2)
g1row<-plot_grid(g1row,legend1,rel_widths = c(4,0.5))
plot(g1row)
dev.off()
png(paste0(result_path,"celltype_umap.png"),width=30, height = 10, units="in", res=100)
g2row<-plot_grid(plotlist=gc, labels =c("i","j","k","l","m","n","o","p"),label_size = 40,nrow = 2)
g2row<-plot_grid(g2row,legend2,rel_widths = c(4,0.5))
plot(g2row)
dev.off()
#png(paste0(result_path,"cell_alignment_umap.png"),width=32, height = 8, units="in", res=100)
#cp<-plot_grid(plotlist = ga, labels = "auto", nrow=2)
#cp2<-plot_grid(cp,plot_grid(legend1,legend2,ncol=1), rel_widths = c(9,0.5))
#dev.off()

# plot metric
g1<-ggplot(data = ari_data, aes(x=group,y=cari))+
	geom_bar(stat="identity",fill="orange")+
	ylab("ARI")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())

# Compute or Extract Silhouette Information
g2<-ggplot(data = sil_data, aes(x=group,y=csil))+
	geom_boxplot()+
	ylab("Silhouette")+
	theme(axis.title.x = element_blank())

# kBET line plot
kbet_data_sub<-kbet_data[kbet_data$neighbor_size==0.5,]
# g5<-ggplot(data= kbet_data_sub,aes(x=cell_type,y = mkb))+
# 	geom_line(aes(color=group))+
# 	ylab("kBET acceptance rate")+
# 	xlab("cell types")+
# 	theme(legend.title = element_blank())
# legend3 <- cowplot::get_legend(g5)
# g5<-g5+theme(legend.position = "none")
#
kbet_tt<-as.data.frame(table(kbet_data_sub[kbet_data_sub$mkb>0.7,]$group))
g3<-ggplot(data = kbet_tt, aes(x=Var1,y=Freq))+
	geom_bar(stat = "identity",fill="orange")+
	ylab("# cell types")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())
png(paste0(result_path,"kbet_bars.png"),width=34,height=2, units="in",res=100)
g3
dev.off()
df=data.frame(mm=mm,aa=c("DESC","LIGER","MNN","scAlign",
												 "Scanorama","scVI","Seurat V3","VIPCCA"))
g4<-ggplot(data = df,aes(x=aa,y=mm))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Mixing Metric")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())

# kBET test metric
png(paste0(result_path,"kbet_celltype.png"),width=28, height = 20, units="in", res=100)
plot_grid(plotlist = gd, labels = "auto", ncol = 1)
dev.off()
g5<-ggplot(data = kbet_data, aes(x=group,y=mkb))+
	geom_boxplot()+
	ylab("kBET (acceptance rate)")+theme_bw(size=30)+
	theme(axis.title.x = element_blank())

library(rlist)
# DE
#rm(list = ls())
source("./code/run/evaluation_function.R")
result_path="./results/pancreas/"
batch_col="dataset"
celltype_col="celltype"
meth_names<-c("MNN","VIPCCA","Scanorama","Seurat V3","scVI")
curr_dirs<-c(
							"mnn/mnn_integrated_data.rds",
							"scxx_cvae/cvae_var2_10_128_16/integrated_data.rds",
						 "scanorama/integrated_data.rds",
						 "seurat/seurat_v3_integrated_data.rds",
						 "scvi/integrated_data.rds"
						 )
br.all<-data.frame()
i=0
for(cdir in curr_dirs){
	i<-i+1
	rfile<-paste0(result_path,cdir)
	#print(strsplit(cdir,split = "/")[[1]][1])
	if(file.exists(rfile))
		print(rfile)
	#browser()
	data.integrated <- readRDS(rfile)
	dim(data.integrated)
	print(length(VariableFeatures(data.integrated)))
	# for each cell type
	# br.all<-foreach(ct=levels(factor(data.integrated@meta.data[[celltype_col]])),
	# 				.export = 'run_dge_analyse',
	# 				.combine = 'rbind')%dopar% run_dge_analyse(data.integrated,
	# 																									 batch1 = "indrop3",
	# 																									 batch2 = "smartseq2",
	# 																									 celltype1 = ct,
	# 																									 celltype_col = celltype_col,
	# 																									 batch_col = batch_col,
	# 																									 method = strsplit(cdir,split = "/")[[1]][1])

	
	for(ct in
			levels(factor(data.integrated@meta.data[[celltype_col]]))
			){
		br<-run_dge_analyse(data.integrated,
												batch1 = "indrop3",
												batch2 = "smartseq2",
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
saveRDS(br.all,file = "./results/pancreas/result_R/sig_genes_all.rds")
br.all<-readRDS("./results/pancreas/result_R/sig_genes_all.rds")
br.all<-br.all[br.all$method!="liger",]
br.all<-br.all[br.all$method!="desc",]
br.all<-br.all[br.all$method!="scalign",]
br.selected <- br.all[br.all$p_val_adj<1e-50,]
sig_genes<-as.data.frame(table(br.selected$method,br.selected$ctype))
sig_genes$Var1<-as.character(sig_genes$Var1)
# c("MNN","Scanorama","scVI","Seurat V3","VIPCCA")
sig_genes$Var2<-as.character(sig_genes$Var2)
sig_genes<-rbind(sig_genes,c("VIPCCA","alpha",0),
								 c("VIPCCA","beta",0),
								 c("VIPCCA","ductal",0),
								 c("VIPCCA","acinar",0),
								 c("VIPCCA","delta",0),
								 c("VIPCCA","activated_stellate",0))
sig_genes$Freq<-as.numeric(sig_genes$Freq)

#source("./code/run/qqplot.R")
# df_qq<-data.frame()
# ct="beta"
# for(md in levels(factor(br.all$method))){
# 	qq<-createData(md,br.all[br.all$ctype==ct & br.all$method==md,"p_val"])
# 	df_qq<-rbind(df_qq,qq)
# }
# p <- ggplot(df_qq, aes(x=expected,y=observed,group=Method))
# p<-p+geom_ribbon(aes(ymin=clower,ymax=cupper),fill = "grey70")+
# 	geom_line(aes(colour=Method))
# +theme(legend.position = c(0.1,0.8))+
# 	theme(legend.title = element_blank())


g6<-ggplot(data = sig_genes,aes(x = Var1, y=Freq,fill=Var2))+
	geom_bar(position="stack", stat="identity")+
	ylab("# significant genes")+
	labs(fill = "cell type")+theme_bw(base_size=30)+
	theme(axis.title.x = element_blank())+
	theme(axis.text.x = element_text(size=20))+
	theme(legend.text=element_text(size=40))+
	theme(legend.position = c(0.55,0.6),legend.title = element_blank())
legend3<-cowplot::get_legend(g6)
g6<-g6+theme(legend.position = "none")

g6
source("./code/run/evaluation_function.R")
br.two<-data.frame()
br2.two<-data.frame()
i<-0
for(cdir in curr_dirs){
	i<-i+1
	rfile<-paste0(result_path,cdir)
	if(file.exists(rfile))
		print(rfile)
	data.integrated <- readRDS(rfile)
	# for each cell type
	
	br<-run_dge_analyse(data.integrated,
											batch1 = "indrop3",
											celltype1 = "alpha",
											celltype2 = "beta",
											celltype_col = celltype_col,
											batch_col = batch_col,
											method = meth_names[i])
	br2<-run_dge_analyse(data.integrated,
											batch1 = "smartseq2",
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
br.selected<- br[br$p_val_adj<1e-50,]
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

br.two.selected <- br.two[br.two$p_val_adj<1e-50,]
br2.two.selected <- br2.two[br2.two$p_val_adj<1e-50,]

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
	theme(legend.position=c(0.2,0.6),legend.title = element_blank(),
				legend.text=element_text(size=40))
legend4<-cowplot::get_legend(g7)
g7<-g7+NoLegend()

g8<-ggplot(data = data.frame(sim=j_index,mds=mds),aes(x = mds, y=sim))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Jaccard Index")+
	theme_bw(base_size = 30)+
	theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20))


g3row<-plot_grid(plotlist = list(g3,g4,g1,g6,legend3,g8,g7,legend4),
								 labels=c("q","r","s","t","","u","v",""),ncol=5,label_size = 40,rel_widths = c(1,1,1,1,0.5))

png(paste0(result_path,"all_in_one.png"),width = 40,height = 40,units = "in",res = 200)
plot_grid(plotlist = list(g1row,g2row,g3row),ncol = 1)
dev.off()


