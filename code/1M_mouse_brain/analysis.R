# cpu and mem usage of vipcca and scanorama
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
# f_vip<-"./results/1M_mouse_brain/cpu_mem_profile_vipcca.txt"
# f_scanorama<-"./results/1M_mouse_brain/cpu_mem_profile_scanorama.txt"
# f_vip_loss<-"./results/1M_mouse_brain/VIPCCA/run-cvae_cvae_2_64_16_logs20200705-013310-tag-loss.csv"
# ann_saunders_url<-"https://storage.googleapis.com/dropviz-downloads/static/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS"

#ann<-read.csv("./results/1M_mouse_brain/VIPCCA/cvae_2_64_16/annotation.csv")
# 
# data_path<-"./data/scanorama_data/data/mouse_brain/"
# bnames<-dir("./data/scanorama_data/data/mouse_brain/dropviz/")
# 
# library(Matrix)
# ann_dropviz<-data.frame()
#for(bname in bnames){
# 	url_outcome<-paste("https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60",bname,".cell_cluster_outcomes.RDS",sep = "")
# 	url<-paste("https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60",bname,".cluster.assign.RDS",sep = "")
# 	url_dge<-paste("https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60",bname,".raw.dge.txt.gz",sep="")
# 	
# 	destfile<-paste("./data/scanorama_data/data/mouse_brain/dropviz/",bname,"/annotation.rds",sep = "")
# 	destfile_outcome<-paste("./data/scanorama_data/data/mouse_brain/dropviz/",bname,"/annotation_outcome.rds",sep = "")
# 	destfile_dge<-paste("./data/scanorama_data/data/mouse_brain/dropviz/",bname,"/data.raw.dge.txt",sep = "")
# 	outfile<-paste("./data/scanorama_data/data/mouse_brain/dropviz/",bname,"/data_subset.h5ad",sep = "")
# 	barcodefile<-paste("./data/scanorama_data/data/mouse_brain/dropviz/",bname,"/barcodes.tsv",sep = "")
# 	# download.file(url_outcome,destfile = destfile_outcome)
# 	# download.file(url,destfile = destfile)
# 	# download.file(url_dge,destfile = destfile_dge)
# 	# print(c(destfile,outfile,destfile_dge))
# 	# print(length(readRDS(destfile)))
# 	# print(dim(readRDS(destfile_outcome)))
	# ann_outcome=readRDS(destfile_outcome)
	# ann_dropviz<-rbind(ann_dropviz,ann_outcome)
# 	
# 	# lines <- grep('CELL_BARCODES',readLines(destfile_dge),value = T)
# 	# rx <- lapply(strsplit(lines,split = "\t"),function(x) x[-1])
# 	# barcodes <- unlist(rx)
# 	# write.table(barcodes,file = barcodefile,row.names = F,col.names = F,quote = F)
# 
# }
# saveRDS(ann_dropviz,file = "./results/1M_mouse_brain/ann_dropviz.rds")
# ann<-readRDS("./data/scanorama_data/data/mouse_brain/dropviz/Cerebellum_ALT/annotation.rds")
# adata<-ReadH5AD("./data/scanorama_data/data/mouse_brain/dropviz/Cerebellum_ALT/data_subset.h5ad")
# bnames<-c(bnames,"spinalCord")
# download.file(ann_saunders_url,destfile = "./data/scanorama_data/data/mouse_brain/ann_saunders.rds")
# ann_saunders<-readRDS("./data/scanorama_data/data/mouse_brain/ann_saunders.rds")
# ann_nuclei<-read.csv("./data/scanorama_data/data/mouse_brain/nuclei/barcodes.txt",
# 										 header = FALSE,sep = '\t',row.names = 1,
# 										 col.names = c("barcode","sample_type","cell_type","spinal_cell_type"))
# saveRDS(ann_nuclei,"./results/1M_mouse_brain/ann_nuclei.rds")


# read barcode subset
set.seed(100)
files<-dir("./data/scanorama_data/data/mouse_brain",pattern="tab.barcodes.txt",recursive=T,full.names = T)
# barcodes_dropviz<-c()
# for(f in files[c(1:9)]){
# 	print(f)
# 	barcodes<-read.csv(f,header = F,stringsAsFactors = F)
# 	barcodes_dropviz<-c(barcodes_dropviz,barcodes)
# }
# barcode_split_seq<-read.csv(files[10],header = F,stringsAsFactors = F)
# barcode_all<-c(barcode_split_seq$V1,unname(unlist(barcodes_dropviz)))
# barcode_dropviz_selected<-sample(unname(unlist(barcodes_dropviz)),
# 																 size = 20000,replace = F)
# barcode_split_seq_selected<-as.character(sample(barcode_split_seq$V1, size = 20000,
# 																	replace = F))
# write.table(barcode_all,file = "./results/1M_mouse_brain/barcode_all.txt",col.names = F,
# 						row.names = F)
# write.table(c(barcode_split_seq_selected,barcodes_dropviz_selected),file = "./results/1M_mouse_brain/barcode_selected.txt",col.names = F,
#						row.names = F)

barcodes_dropviz_selected<-c()
for(f in files[c(1:9)]){
	print(f)
	barcodes<-read.csv(f,header = F,stringsAsFactors = F)
	barcodes_selected<-sample(barcodes$V1,size = dim(barcodes)[1]/10)
	barcodes_dropviz_selected<-c(barcodes_dropviz_selected,barcodes_selected)
}
barcode_split_seq<-read.csv(files[10],header = F,stringsAsFactors = F)
barcode_split_seq_selected<-as.character(sample(barcode_split_seq$V1, 
																								size = dim(barcode_split_seq)[1]/10,
																								replace = F))
barcode_all <- read.csv("./results/1M_mouse_brain/barcode_all.txt",header = F,stringsAsFactors = F)$V1
# write.table(barcode_all,file = "./results/1M_mouse_brain/barcode_all.txt",col.names = F,
# 						row.names = F)
write.table(c(barcode_split_seq_selected,barcodes_dropviz_selected),file = "./results/1M_mouse_brain/barcode_selected_2.txt",col.names = F,
						row.names = F)

# # subset of scanorama output
# adata_scanorama<-ReadH5AD("./results/1M_mouse_brain/Scanorama/output_unq.h5ad")
# adata_scanorama<-RenameCells(adata_scanorama, new.names = barcode_all)
# adata_scanorama_subset<-subset(adata_scanorama,cells=c(barcode_split_seq_selected,
# 																											 barcodes_dropviz_selected))
# saveRDS(adata_scanorama_subset,
# 				file = "./results/1M_mouse_brain/Scanorama/output_unq_subset_2.rds")
# rm(adata_scanorama,adata_scanorama_subset)
# 
# # subset of vipcca output
# f_vip_data<-"./results/1M_mouse_brain/VIPCCA/cvae_2_64_16/output.h5ad"
# adata_vipcca<-ReadH5AD(f_vip_data)
# adata_vipcca_subset<-subset(adata_vipcca,cells=c(barcode_split_seq_selected,
# 																											 barcodes_dropviz_selected))
# saveRDS(adata_vipcca_subset,
# 				file = "./results/1M_mouse_brain/VIPCCA/output_unq_subset_2.rds")
# rm(adata_vipcca,adata_vipcca_subset)

# # subset of harmony output
f_harmony_data<-"./results/1M_mouse_brain/harmony/integrated_data.rds"
adata_harmony<-readRDS(f_harmony_data)
adata_harmony_subset<-subset(adata_harmony,cells=c(barcode_split_seq_selected,
																											 barcodes_dropviz_selected))
saveRDS(adata_harmony_subset,
				file = "./results/1M_mouse_brain/harmony/output_unq_subset_2.rds")
rm(adata_harmony,adata_harmony_subset)

# 
# # subset of others
# library(Seurat)
# adata_desc<-ReadH5AD("./results/1M_mouse_brain/desc/output.h5ad")
# tail(colnames(adata_desc))
# 
# adata_desc<-RenameCells(adata_desc, new.names = barcode_all)
# adata_desc_subset<-subset(adata_desc,cells=c(barcode_split_seq_selected,
# 																						barcodes_dropviz_selected))
# saveRDS(adata_desc_subset,
# 				file="./results/1M_mouse_brain/desc/output_unq_subset_2.rds")
# rm(adata_desc,adata_desc_subset)
# # adata<-adata_desc_subset
# #adata_scvi<-ReadH5AD("./results/1M_mouse_brain/scVI/output.h5ad")
# 
# barcode_all<-read.csv("./results/1M_mouse_brain/barcode_all.txt",header = F,stringsAsFactors = F)$V1
# barcode_selected<-read.csv("./results/1M_mouse_brain/barcode_selected_2.txt",header = F,stringsAsFactors = F)$V1
# intersect(barcode_all,barcode_selected)
###
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
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


analyse_on_mouse_brain<-function(adata=NULL, method="VIPCCA",reduction="scxx",
																 result_path="./"){
	# browser()
	result_path<-paste0("./results/1M_mouse_brain/",method,"/", sep = "")
	if(!file.exists(paste0(result_path,method,"_umap_2.rds",sep = ""))){
		adata$.batch<-as.factor(adata$.batch)
		levels(adata$.batch)<-c("nuclei","Cerebellum_ALT",
														"Cortex_noRep5_FRONTALonly","Cortex_noRep5_POSTERIORonly",
														"EntoPeduncular","GlobusPallidus",
														"Hippocampus","Striatum",
														"SubstantiaNigra","Thalamus")
		adata@meta.data["tech"]="drop-seq"
		adata$tech[adata$.batch=="nuclei"]="SPLiT-seq"
		adata <- FindNeighbors(adata, reduction=reduction,dims = 1:16)
		#adata<-FindClusters(adata,resolution = 0.3)
		adata<-RunUMAP(adata,reduction=reduction,dims=1:16)
		saveRDS(adata,paste0(result_path,method,"_umap_2.rds",sep=""))
	}else{
		adata<-readRDS(paste0(result_path,method,"_umap.rds",sep=""))
	}
	
	Idents(adata)<-".batch"
	adata_splitseq<-subset(adata, idents = "nuclei")
	adata_splitseq<-AddMetaData(adata_splitseq,metadata = ann_nuclei)
	adata_dropviz<-adata[,Idents(adata)!="nuclei"]
	adata_dropviz<-AddMetaData(adata_dropviz, metadata = ann_dropviz)
	g1<-DimPlot(adata,reduction = "umap",group.by = "tech")+ggtitle(method)+
		theme_bw(base_size=30)+NoLegend()+theme(axis.text.x = element_text(size=15))+
		theme(axis.title.x = element_blank())
	g2<-DimPlot(adata_dropviz,reduction = "umap",group.by = "cluster")+
		ggtitle(paste(method,"drop-seq"))+theme_bw(base_size=30)+
		NoLegend()+theme(axis.text.x = element_text(size=15))+
		theme(axis.title.x = element_blank())
	g3<-DimPlot(adata_splitseq,reduction = "umap",group.by = "cell_type")+
		ggtitle(paste(method,"SPLiT-seq"))+theme_bw(base_size=30)+
		NoLegend()+theme(axis.text.x = element_text(size=15))+
		theme(axis.title.x = element_blank())
	
	# ARI for seurat clusters
	# browser()
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
	kfile=paste0("kbet_",method,"_2.rds",sep="")
	celltype_col="is_selected"
	batch_col="tech"
	if(!file.exists(paste0(result_path,kfile,sep=""))){
		require('FNN')
		com=run_kbet_for_each_celltype2(adata,
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
	g4<-com[[2]]+ylim(c(0,0.4))+theme_bw(base_size=30)+theme(axis.text.x = element_text(size=15))+
		theme(axis.title.x = element_blank())+ylab("kBET")
	
	#browser()
	max.k=300
	mm <- max.k-MixingMetric(adata, grouping.var=batch_col, 
													 reduction=reduction, dims=1:16, max.k=max.k)

	return(list(g1,g2,g3,g4,median(mm),df_ari))
}

adata_vipcca<-readRDS("./results/1M_mouse_brain/VIPCCA/output_unq_subset_2.rds")
adata_scanorama<-readRDS("./results/1M_mouse_brain/Scanorama/output_unq_subset_2.rds")
adata_scvi <- ReadH5AD("./results/1M_mouse_brain/scVI/output_sub_2.h5ad")
adata_desc <-  readRDS("./results/1M_mouse_brain/desc/output_unq_subset_2.rds")

glist_vipcca<-analyse_on_mouse_brain(adata = adata_vipcca, reduction = "scxx", method="VIPCCA")
glist_scanorama<-analyse_on_mouse_brain(adata = adata_scanorama, method = "Scanorama", reduction = "scanorama")
glist_scvi <- analyse_on_mouse_brain(adata = adata_scvi,	method="scVI", reduction = "scvi")
glist_desc <- analyse_on_mouse_brain(adata = adata_desc,	method="DESC", reduction = "Embeded_z0.3")

library(rlist)
df<-data.frame(mm=c(glist_vipcca[[5]],glist_scanorama[[5]],glist_scvi[[5]],glist_desc[[5]]),
					 method=c("VIPCCA","Scanorama","scVI","DESC"))
p1<-ggplot(data = df,aes(x=method,y=mm))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Mixing Metric")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())

df_ari<-do.call(rbind, list(glist_vipcca[[6]],glist_scanorama[[6]],glist_scvi[[6]],glist_desc[[6]]))

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
f_scvi <- paste("./results/1M_mouse_brain/cpu_mem_profile_scvi_new.txt", sep = "")
x3<-read.delim(f_scvi, header = F, sep = ",", comment.char = "#",col.names = cnames,stringsAsFactors = T)
x3$method="scVI"
x3$ctime=(0:(dim(x3)[1]-1))*1.0/3600
x3$mem_In_G<-x3$mem*2.56
f_desc <- paste("./results/1M_mouse_brain/cpu_mem_profile_desc_new.txt", sep = "")
x4<-read.delim(f_desc, header = F, sep = ",", comment.char = "#",col.names = cnames,stringsAsFactors = T)
x4$method="DESC"
x4$ctime=(0:(dim(x4)[1]-1))*1.0/3600
x4$mem_In_G<-x4$mem*2.56

df_scale <- do.call(rbind, list(x1,x2,x3,x4))

df_cpu_mean<-aggregate(cpu~method,df_scale,FUN = mean)
df_cpu_peak<-aggregate(cpu~method,df_scale,FUN = max)
df_mem_mean<-aggregate(mem_In_G~method,df_scale,FUN = mean)
df_mem_peak<-aggregate(mem_In_G~method,df_scale,FUN = max)
df_time_used<-aggregate(ctime~method,df_scale,FUN = max)

p2<-ggplot(data = df_cpu_mean,aes(x=method,y=cpu))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Average cpu (%)")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())
p3<-ggplot(data = df_cpu_peak,aes(x=method,y=cpu))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Peak cpu (%)")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())
p4<-ggplot(data = df_mem_mean,aes(x=method,y=mem_In_G))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Average mem (G)")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())
p5<-ggplot(data = df_mem_peak,aes(x=method,y=mem_In_G))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Peak mem (G)")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())
p6<-ggplot(data = df_time_used,aes(x=method,y=ctime))+
	geom_bar(stat="identity",fill="orange")+
	ylab("Time elapsed (H)")+theme_bw(base_size=30)+
	theme(axis.text.x = element_text(size=15))+
	theme(axis.title.x = element_blank())


glist_all<-c(glist_desc[1:4],glist_scanorama[1:4],glist_scvi[1:4],glist_vipcca[1:4])
glist_all<-glist_all[c(seq(1,16,4),seq(2,16,4),seq(3,16,4),seq(4,16,4))]
glist_all<-list.append(glist_all,p1,p2,p3,p4,p5,p6)
png("./results/1M_mouse_brain/plot_all.png",width = 24, height = 24,units = "in", res = 300)
plot_grid(plotlist = glist_all,nrow = 6,labels = "auto",label_size = 30)
dev.off()

rm(list = ls())
library(rlist)
methods <- c("DESC", "Scanorama", "scVI", "VIPCCA")
sp_list <- list()
for(method in methods){
	adata<-readRDS(paste0("./results/1M_mouse_brain/",method,"/ ",method," _umap.rds",sep=""))
	levels(adata@meta.data$.batch)[1]="CNS (SPLiT-seq)"
	sp <- DimPlot(adata,reduction = "umap",group.by = ".batch")+ggtitle(method)+
		theme_bw(base_size=30)+theme(axis.text.x = element_text(size=15))+
		theme(axis.title.x = element_blank())
	legend_region<-cowplot::get_legend(sp)
	sp <- sp + NoLegend()
	sp_list<-list.append(sp_list,sp)
}
png("./results/1M_mouse_brain/sp.png",width = 16, height = 8, units = "in", res = 300)
plot_grid(plot_grid(plotlist = sp_list),legend_region, rel_widths = c(3,2))
dev.off()

## markers ...
#	mouse_brain_genes = [
	#     'Gja1', 'Flt1', 'Gabra6', 'Syt1', 'Gabrb2', 'Gabra1',
	#     'Meg3', 'Mbp', 'Rgs5',
	# ]


# library(ggplot2)
# library(cowplot)
# p1<-ggplot(x,aes(x=ctime,y=mem,colour=method))+
# 	geom_line()+xlab("hours")+ylab("mem (GB)")+
# 	theme(legend.position = c(0.2, 0.8))
# p2<-ggplot(x,aes(x=ctime,y=cpu,colour=method))+
# 	geom_line()+xlab("hours")+ylab("cpu (%)")+
# 	theme(legend.position = c(0.8, 0.8))
# p3<-ggplot(x3,aes(x=ctime,y=Value,colour=method))+
# 	#xlim(0,max(x$ctime))+
# 	geom_line()+xlab("hours")+ylab("loss")+
# 	theme(legend.position = c(0.8,0.8))
# 
# g3<-plot_grid(g1,g2,nrow = 2,labels = c("a","b"))
# g4<-plot_grid(p1,p2,p3,nrow = 3,labels = c("c","d","e"))
# g5<-plot_grid(g3,g4,nrow=1)
# png("./results/1M_mouse_brain/cpu_mem_loss_plot.png",width = 18, height = 18, units = "in", res = 300)
# g5
# dev.off()
