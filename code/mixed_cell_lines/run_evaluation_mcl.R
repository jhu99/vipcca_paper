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
result_path="./results/mixed_cell_lines/real/"
celltype_col="celltype"
batch_col="X_batch"
cl <- makeCluster(20)
registerDoParallel(cl)
clusterExport(cl, varlist=c("kbet_for_one_celltype", "get.knn", "kBET"), env=environment())

res_list <- list()
ari<-c()
gari<-c()
rr<-c()
lrr<-c()
gd<-list()
curr_path=paste0(result_path,"vipcca/")
if(!file.exists(paste0(curr_path,"kbet_cvae.png"))){
	methods<-dir(curr_path,pattern = "cvae")
	for(dir_cvae in methods){
		print(dir_cvae)
		h5ad_to_seurat(input_file = "output.h5ad", result_path = paste0(curr_path,dir_cvae,"/"),
									 annotation_file="./data/mixed_cell_lines/annotation.txt")
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
saveRDS(list(ari_data_cvae,kbet_data_cvae),file = "./results/mixed_cell_lines/real/result_R/cvae_stat.rds")
ax1<-readRDS("./results/mixed_cell_lines/real/result_R/cvae_stat.rds")
ad1<-aggregate(cari~group,ax1[[1]],max)
md1<-aggregate(mkb~group,ax1[[2]],median)
md1[md1$mkb>0.9,]
ad1[rownames(md1[md1$mkb>0.9,]),]

#cvae
source("./code/run/evaluation_function.R")
curr_path=paste0(result_path,"vipcca/cvae_2_64_14/")
vipcca_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="VIPCCA",algorithm="cover_tree",limmax=20,
																			 batch_col=batch_col,celltype_col=celltype_col,reduction="scxx",dims=1:16,min.dist=0.01,fmt="rds")
curr_path=paste0(result_path,"liger/")
liger_res <- calculate_alignment_metric(curr_path,
																				"liger_integrated_data.rds",method="LIGER",limmax=20,
																				batch_col=batch_col,celltype_col=celltype_col,reduction="iNMF",dims=1:16,fmt="rds")

curr_path=paste0(result_path,"seurat/")
seurat_res <- calculate_alignment_metric(curr_path,"seurat_v3_integrated_data.rds",method="Seurat V3",limmax=20,
																				 batch_col=batch_col,celltype_col=celltype_col,reduction="pca",dims=1:16,fmt="rds")

curr_path=paste0(result_path,"scanorama/")
scanorama_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="Scanorama",limmax=20,
																						batch_col=batch_col,celltype_col=celltype_col,reduction="scanorama",dims=1:16,fmt="rds")


curr_path=paste0(result_path,"desc/")
desc_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="DESC",limmax=20,
																						batch_col=batch_col,celltype_col=celltype_col,reduction="Embeded_z0.3",dims=1:16,fmt="rds")

# curr_path=paste0(result_path,"scvi/")
# scvi_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="scVI",
# 																			 batch_col=batch_col,celltype_col=celltype_col,reduction="scvi",dims=1:16,fmt="rds" )

curr_path=paste0(result_path,"scalign/")
scalign_res <- calculate_alignment_metric(curr_path,"scalign_seuratobj.rds",method="ScAlign",limmax=20,
																			 batch_col=batch_col,celltype_col=celltype_col,reduction="ALIGNED.MultiCCA",dims=1:16,fmt="rds" )

curr_path=paste0(result_path,"mnn/")
mnn_res <- calculate_alignment_metric(curr_path,"mnn_integrated_data.rds",method="MNN",limmax=20,
																			 batch_col=batch_col,celltype_col=celltype_col,reduction="pca",dims=1:16,fmt="rds" )

curr_path=paste0(result_path,"harmony/")
harmony_res <- calculate_alignment_metric(curr_path,"integrated_data.rds",method="Harmony",limmax=20,
																			batch_col=batch_col,celltype_col=celltype_col,reduction="harmony",dims=1:16,fmt="rds" )


datasets.integrated<-readRDS(paste0(result_path, "/vipcca/cvae_2_64_14/integrated_processed.rds"))
datasets.integrated$X_batch <- as.factor(datasets.integrated$X_batch)
levels(datasets.integrated$X_batch) <- c("Dataset1","Dataset2","Dataset3")
datasets.integrated$celltype <- as.factor(datasets.integrated$celltype)
levels(datasets.integrated$celltype)<-c("293t","Jurkat")

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(3)
p1<-DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col) +
	scale_color_manual(values=mycolors)+
	theme(legend.text=element_text(size=70,face = "bold"),legend.position=c(0.4,0.5))+
	guides(colour = guide_legend(override.aes = list(size=15),nrow=1))
legend1 <- cowplot::get_legend(p1)
p2<-DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col) +
	theme(legend.text=element_text(size=70,face = "bold"),legend.position=c(0.4,0.5))+
	guides(colour = guide_legend(override.aes = list(size=15),nrow=1))
legend2 <- cowplot::get_legend(p2)

library(rlist)
alg_names<-c("DESC","Harmony", "LIGER","MNN","ScAlign","Scanorama","Seurat V3","VIPCCA")

ga<-list()
gb<-list()
gc<-list()
gd<-list()
ge<-list()
ari_data<-c()
sil<-c()
mm<-c()
tt<-c()
rr<-c()
lrr<-c()
df_kbet<-list()
i<- 0
for(res in list(desc_res,harmony_res,liger_res,mnn_res,scalign_res,
								scanorama_res,seurat_res,vipcca_res)){
	#ga<-list.append(ga,res[[4]])
	#
  i <- i+1
	mm <- c(mm,res[[1]])
	sil <- c(sil,res[[2]])
	ari_data <- rbind(ari_data,res[[3]]) 
	tt<- c(tt,length(res[[2]]))
	rr<-c(rr,res[[6]][[1]][,3])
	lrr<-c(lrr,length(res[[6]][[1]][,3]))
	gb<-list.append(gb,res[[4]])
	gc<-list.append(gc,res[[5]])
	gd<-list.append(gd,res[[6]][[2]])
	ge<-list.append(ge,res[[7]][[2]]+ylim(0,1)+ggtitle(alg_names[i]))
	df_kbet<-list.append(df_kbet,res[[7]][[1]])
}
df_kbet<-do.call("rbind",df_kbet)
df_kbet[,"methods"]=rep(alg_names,each=10)


# ari_data$group <- rep(c("DESC","Harmony","LIGER","MNN","scAlign",
# 												"Scanorama","Seurat V3","VIPCCA"),
# 											each=10)
ari_data$group <- rep(c("DESC","Harmony","LIGER","MNN","ScAlign",
                        "Scanorama","Seurat V3","VIPCCA"),each=10)
ari_data$ari_batch<-1-ari_data$ari_batch
ss<-seq(from=1, to = 80,by = 10)
ari_data<-ari_data[ss,]
sil_data<- data.frame(csil=sil,
											group=rep(alg_names,
																				 tt))
kbet_data<- data.frame(mkb=rr,group=rep(alg_names,lrr),
											 neighbor_size=rep(seq(5,25,by = 5),2*4))

png(paste0(result_path,"dataset_umap.png"),width=50, height = 20, units="in", res=100)
g1row<-plot_grid(plotlist=gb, nrow = 2)
# g1row<-g1row+draw_label("UMAP-1",x = 0.4,y=0,hjust = 0,size = 70,fontface = "bold")+
# 	draw_label("UMAP-2",	x = 0,y=0.4,hjust = 0,size = 70,angle = 90,fontface = "bold")+
# 	theme(plot.margin = margin(0, 0, 100, 100,unit = "pt"))
g1row<-plot_grid(g1row,legend1,rel_heights  = c(6,1.0),ncol = 1)
plot(g1row)
dev.off()

png(paste0(result_path,"celltype_umap.png"),width=50, height = 20, units="in", res=100)
g2row<-plot_grid(plotlist=gc, nrow = 2)
# g2row<-g2row+draw_label("UMAP-1",x = 0.4,y=0,hjust = 0,size = 70,fontface = "bold")+
# 	draw_label("UMAP-2",	x = 0,y=0.4,hjust = 0,size = 70,angle = 90,fontface = "bold")+
# 	theme(plot.margin = margin(0, 0, 100, 100,unit = "pt"))
g2row<-plot_grid(g2row,legend2,rel_heights  = c(6,1.0),ncol = 1)
plot(g2row)
dev.off()

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
				legend.text = element_text(size = 70,face="bold"),legend.key.size = unit(4,"cm"),legend.key.height = unit(2,"cm"),
				axis.text = element_text(size = 60,face = "bold"),
				axis.title = element_text(size = 70,face = "bold"),legend.position = c(0.5,0.3),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
	guides(colour = guide_legend(override.aes = list(size=15),nrow = 1))
legend3<-cowplot::get_legend(g5)
g5<-g5+NoLegend()
g3row<-plot_grid(plotlist = list(g4,g5),
								 labels=c("c","d"),nrow=1,label_size = 100,rel_widths = c(2,2),
								 label_x = 0, label_y = 1,
								 hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,3),"cm"))
png(paste0(result_path,"mm_kbet.png"),width=50, height = 15, units="in", res=100)
plot(g3row)
dev.off()


# plot metric
# g1<-ggplot(data = ari_data, aes(x=ari_celltype,y=ari_batch))+
# 	geom_point(aes(group=group,color=group),size=10)+xlab("ARI celltype")+
# 	ylab("1-ARI batch")+theme_bw(base_size=30)+
# 	theme(axis.text.x = element_text(size=15),legend.position = "top",legend.title = element_blank())
g1<-ggplot(data = ari_data, aes(y=ari_celltype,x=group))+
	geom_bar(stat="identity",fill=colors)+
	ylab("ARI")+
	theme(axis.title.x = element_blank(),
				axis.text.y  = element_text(size = 60,face = "bold"),
				axis.text.x  = element_blank(),
				axis.title.y = element_text(size = 70,face="bold"),
				legend.title = element_blank(),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())


# kBET line plot
# kbet_data_sub<-kbet_data[kbet_data$neighbor_size==0.5,]
# kbet_tt<-as.data.frame(table(kbet_data_sub[kbet_data_sub$mkb>0.7,]$group))
# g3<-ggplot(data = kbet_tt, aes(x=Var1,y=Freq))+
# 	geom_bar(stat = "identity",fill="orange")+
# 	ylab("# cell types")+theme_bw(base_size=30)+
# 	theme(axis.text.x = element_text(size=15))+
# 	theme(axis.title.x = element_blank())



png(paste0(result_path,"kbet_celltype.png"),width=8, height = 16, units="in", res=100)
plot_grid(plotlist= gd, labels = "auto", ncol = 1)
dev.off()

png(paste0(result_path,"kbet_celltype2.png"),width=20, height = 35, units="in", res=100)
plot_grid(plotlist= ge, labels = "auto", ncol = 1,label_size = 40)
dev.off()

levels(df_kbet$celltype)<-c("293t","Jurkat")
png(paste0(result_path,"kbet_celltype3.png"),width=15, height = 7, units="in", res=100)
pkbet3<-ggplot(df_kbet,aes(x=nb,y=ymedian,group=methods))+geom_line(aes(color=methods),size=5)+
	scale_color_manual(values=colors)+
	theme(text = element_text(size=40,face = "bold"),panel.grid.major = element_blank(),legend.position = "bottom",
				panel.grid.minor = element_blank(),panel.background = element_blank(),
				legend.title = element_blank(),legend.key.size = unit(2,"cm"),
				axis.line = element_line(colour = "black"))+labs(x="% Sample size",y="kBET")+
	facet_wrap(~celltype,ncol = 2)
dev.off()
# kbet_data_sub = kbet_data[kbet_data$neighbor_size==0.6,]
# g5<-ggplot(data= kbet_data_sub,aes(x=cell_type,y = mkb,fill=group))+
# 	geom_bar(position="dodge",stat="identity")+
# 	ylab("kBET acceptance rate")+
# 	xlab("cell types")+
# 	theme(legend.title = element_blank(),
# 				legend.position = c(0.9,0.9))
# legend3 <- cowplot::get_legend(g5)
# g5<-g5+theme(legend.position = "none")

source("./code/run/evaluation_function.R")
result_path="./results/mixed_cell_lines/real/"
batch_col="X_batch"
celltype_col="celltype"
meth_names<-c("VIPCCA","Seurat V3","Scanorama","MNN")
vipcca_adata<-readRDS(paste0(result_path,"vipcca/cvae_2_64_14/integrated_data.rds"))
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
  i <- i+1
	# for each cell type
	# foreach(ct=levels(factor(data.integrated@meta.data[[celltype_col]])),
	# 				.export = 'run_dge_analyse',
	# 				.combine = 'rbind')%dopar% run_dge_analyse(data.integrated,
	# 																									 batch1 = "indrop3",
	# 																									 batch2 = "smartseq2",
	# 																									 celltype1 = ct,
	# 																									 celltype_col = celltype_col,
	# 																									 batch_col = batch_col,
	# 																									 method = strsplit(cdir,split = "/")[[1]][1])
	# 
	for(ct in levels(factor(data.integrated@meta.data[[celltype_col]]))
	){
	  browser()
		br<-run_dge_analyse(data.integrated,
												batch1 = "mixed",
												batch2 = ct,
												celltype1 = ct,
												celltype_col = celltype_col,batch_col = batch_col,method = meth_names[i])
		if(!is.null(br)){
			br.all<-rbind(br.all,br)
		}
		else{
			print("Cells are less than 50 in at least one group.")
		}
		print(ct)
	}
}
saveRDS(br.all,"./results/mixed_cell_lines/real/result_R/sig_genes_all.rds")
br.all<-readRDS("./results/mixed_cell_lines/real/result_R/sig_genes_all.rds")
br.all$log_p_val<- -log10(br.all$p_val+1e-300)
br.all[which(br.all$ctype=="jurkat"),"ctype"]="Jurkat"
g6<-ggplot(br.all, aes(y=log_p_val,x=ctype,fill=method))+
	geom_boxplot()+
	scale_fill_manual(values = colors[c(4,6,7,8)])+
	labs(y="-log (p-value)")+
	theme(legend.title = element_blank(),
				legend.text = element_text(size = 70,face = "bold"),
				axis.text = element_text(size = 60,face = "bold"),
				axis.title.x = element_blank(),
				axis.title = element_text(size = 70,face = "bold"),legend.position = c(0.3,0.6),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
	theme(legend.key.size = unit(2,"cm"))

legend4 <- cowplot::get_legend(g6)
g6<-g6+NoLegend()


g4row<-plot_grid(plotlist = list(g1,g6),
								 labels=c("e","f",""),ncol=2,label_size = 100,rel_widths = c(2,2),
								 label_x = 0, label_y = 1,
								 hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,3),"cm"))
png(paste0(result_path,"ari_pvalue.png"),width=50, height = 15, units="in", res=100)
plot(g4row)
dev.off()
png(paste0(result_path,"all_in_one.png"),width = 50,height = 70,units = "in",res = 100)
plot_grid(plotlist = list(g1row,g2row,g3row,g4row,legend3),ncol = 1,rel_heights =  c(2,2,1.5,1.5,0.2),labels = c("a","b",""),label_size = 100)
dev.off()


aggregate(log_p_val~method,data=br.all,median)

# br.selected <- br.all[br.all$p_val_adj<5e-2,]
# sig_genes<-as.data.frame(table(br.selected$method,br.selected$ctype))
# sig_genes$Var1<-as.character(sig_genes$Var1)
# 
# # sig_genes<-rbind(sig_genes,data.frame(Var1=c("VIPCCA","VIPCCA"),Var2=c("293t","jurkat"),Freq=c(0,0)))
# 
# g6<-ggplot(data=sig_genes,aes(x=Var2,y=Freq))+
#   geom_point(aes(color=Var1),size=8)+theme_bw(base_size=30)+
#   labs(y="# significant genes")+
#   theme(axis.title.x = element_blank())+
#   theme(axis.text.x = element_text(size=20))+
#   theme(legend.title =element_blank())+
#   theme(legend.position = c(0.3,0.6))

legend3<-cowplot::get_legend(g5)
g5<-g5+theme(legend.position = "none")


br.two<-data.frame()
br2.two<-data.frame()
i<-0
source("./code/run/evaluation_function.R")
for(data.integrated in c(vipcca_adata,
                         seurat_adata,
                         scanorma_adata,
                         mnn_adata
)){
	i<-i+1
	# for each cell type
	
	br<-run_dge_analyse(data.integrated,
											batch1 = "mixed",
											celltype1 = "293t",
											celltype2 = "jurkat",
											celltype_col = celltype_col,
											batch_col = batch_col,
											method = meth_names[i])
	br2<-run_dge_analyse(data.integrated,
											 batch1 = "293t",
											 batch2 = "mixed",
											 celltype1 = "293t",
											 celltype2 = "jurkat",
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
saveRDS(br.two,"./results/mixed_cell_lines/real/result_R/br_two_2.rds")
saveRDS(br2.two,"./results/mixed_cell_lines/real/result_R/br2_two_2.rds")
br.two<-readRDS("./results/mixed_cell_lines/real/result_R/br_two_2.rds")
br2.two<-readRDS("./results/mixed_cell_lines/real/result_R/br2_two_2.rds")
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
	theme(legend.title = element_blank(),legend.position = c(0.2,0.9),
				legend.text=element_text(size=30),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank())

#legend4<-cowplot::get_legend(g7)
g8<-ggplot(data = data.frame(sim=j_index,mds=mds),aes(x = mds, y=sim))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Jaccard Index (significant genes)")+
	theme_bw(base_size = 30)+
	theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20))

g3row<-plot_grid(plotlist = list(g4,g5,legend3),
								 labels=c("c","d",""),nrow=1,label_size = 60,rel_widths = c(2.5,1.5,0.5))
g4row<-plot_grid(plotlist = list(g1,g6,legend4),
								 labels=c("e","f",""),nrow=1,label_size = 60,rel_widths = c(2.5,1.5,0.5))


png(paste0(result_path,"plot.all.png"),width = 40, height = 40, units = "in", res=100)
plot_grid(plotlist = list(g1row,g2row,g3row,g4row),rel_heights = c(2,2,1,1),ncol = 1,labels = c("a","b",""),label_size = 60)
dev.off()
##

br.two<-readRDS("./results/mixed_cell_lines/real/result_R/br_two_2.rds")
br2.two<-readRDS("./results/mixed_cell_lines/real/result_R/br2_two_2.rds")
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
	theme(axis.title.x = element_blank(),axis.text = element_text(size=40,face="bold"),
				axis.title.y = element_text(size=40,face="bold"),
				panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())

png(paste0(result_path,"sp_1.png"),width = 20,height = 20,units = "in",res = 100)
plot_grid(pkbet3,g9,labels = "auto",ncol = 1,label_size = 40,label_x = 0, label_y = 1,
					hjust = 0.5, vjust = 0.5)+theme(plot.margin=unit(c(1.5,0,0,3),"cm"))
dev.off()
###########
rm(list = ls())
br.two<-readRDS("./results/mixed_cell_lines/real/result_R/br_two_2.rds")
br2.two<-readRDS("./results/mixed_cell_lines/real/result_R/br2_two_2.rds")
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
	# 								 qvalueCutoff=0.05))
	# kk3 <- as.data.frame(enrichKEGG(geneset3, organism="hsa", pAdjustMethod="BH", 
	# 									qvalueCutoff=0.05))
	kk <- as.data.frame(enrichKEGG(geneset, organism="hsa", pAdjustMethod="BH", 
	                                qvalueCutoff=0.05))
	
	
	
	# kk2$report <- paste0(mt, "-s2")
	# kk3$report <- paste0(mt, "-s3")
	kk$report <- mt
	
	#kk <- rbind(kk2,kk3)
	#df_kk<-rbind(df_kk,rbind(kk2,kk3))
	df_kk <- rbind(df_kk,kk)
	
	# iset <- intersect(kk2$ID,kk3$ID)
	# uset <- union(kk2$ID,kk3$ID)
	# df_kk_sp1<-rbind(df_kk_sp1,list(length(kk2$ID),length(kk3$ID),(1.0*length(iset))/length(uset)))
	
}

df_kk$tqvalue <- -log10(df_kk$qvalue)
st_kk<-as.data.frame(table(df_kk$Description))
df_kk$highlight <- "conserved"
df_kk[df_kk$Description %in% as.character(st_kk[st_kk$Freq<2,"Var1"]),"highlight"]="interest"


# df_kk$highlight <- "no"
# df_kk[df_kk$Description=="Apoptosis","highlight"]="yes"
# df_kk[df_kk$Description=="FoxO signaling pathway","highlight"]="yes"
# df_kk[df_kk$Description=="Mismatch repair","highlight"]="yes"
# df_kk[df_kk$Description=="Wnt signaling pathway","highlight"]="yes"
# df_kk[df_kk$Description=="Transcriptional misregulation in cancer","highlight"]="yes"

g9 <- ggplot(df_kk, aes(y=reorder(Description,tqvalue),x=tqvalue,fill=highlight))+
	geom_bar(stat = "identity")+
	scale_x_continuous(trans = "log10")+labs(x="-log(qvalue)")+
	theme_bw(base_size = 60)+
	theme(axis.title.y = element_blank(),legend.position = "none")+
	facet_wrap(facets = .~report,nrow = 1)+
  ggtitle("Enriched KEGG pathway between 293T and Jurkat")

result_path="./results/mixed_cell_lines/real/"
png(paste0(result_path,"sp_2.png"),width = 55,height = 30,units = "in",res = 100)
g9
dev.off()

colnames(df_kk_sp1)<-c("Scenario2","Scenario3","Jaccard")
df_kk_sp1$method <-  c("VIPCCA","MNN","Scanorama","Seurat V3")


g10 <-ggplot(df_kk_sp1, aes(x=method,y=Jaccard))+
  geom_bar(stat = "identity",fill="orange") +
  ylab("Jaccard Index (enriched pathways)")+
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 30),
        legend.position = c(0.2,0.8),
        axis.title.x = element_blank(),
        legend.text  = element_text(size = 20),legend.title = element_text(size = 30))


png(paste0(result_path,"sp_1.png"),width = 30,height = 10,units = "in",res = 100)
plot_grid(plotlist = list(g7,g8),nrow=1,labels = c("a","b"),label_size = 40)
dev.off()
