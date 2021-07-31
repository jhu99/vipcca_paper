library(Seurat)
library(ggplot2)
library(mclust)
library(cluster)
library(cowplot)
#library(kBET)
library(foreach)
library(wesanderson)
#library(monocle3)
#library(EnhancedVolcano)

calculate_alignment_metric <- function(result_path,
																			 input_file,
																			 method="mnn",
																			 celltype_col="celltype",
																			 batch_col="dataset",
																			 reduction="pca",
																			 dims=1:30,
																			 fmt="h5ad",
																			 min.dist=0.3,
																			 downsample=FALSE,
																			 run_dge=FALSE,
																			 algorithm="cover_tree",
																			 res=0.3,
																			 limmax=15,
																			 integrated_file="integrated_processed.rds"
																			){
	# browser()
	if(!file.exists(paste0(result_path,integrated_file))){
		if(fmt=="h5ad"){
			datasets.integrated=ReadH5AD(paste0(result_path,input_file),assay = "integrated")
			datasets.integrated@assays$integrated@var.features=datasets.integrated@assays$integrated@data@Dimnames[[1]]
		}
		else{
			datasets.integrated <- readRDS(paste0(result_path,input_file))
		}
		# if(method=="scalign"){
		# 	datasets.integrated = as.Seurat(datasets.integrated, counts="counts", scale.data="scale.data")
		# }
		datasets.integrated <- RunUMAP(datasets.integrated, reduction = reduction, dims = dims, min.dist = min.dist, n_threads=10, scale=TRUE)
		datasets.integrated <- FindNeighbors(datasets.integrated, dims = dims, reduction = reduction)
		saveRDS(datasets.integrated, paste0(result_path,integrated_file))
	}else{
		datasets.integrated <- readRDS(paste0(result_path,integrated_file))
	}
	
	#browser()
	celltypes <- as.factor(datasets.integrated@meta.data[[celltype_col]])
	batches <- datasets.integrated@meta.data[[batch_col]]
	tc<-data.frame(batches=batches, celltypes=celltypes)
	
	batch_type = unique(batches)
	gbatch<-ggplot(tc,aes(batches))+geom_bar(aes(fill=celltypes))+theme(legend.position = "right")
	gcelltype<-ggplot(tc,aes(celltypes))+geom_bar(aes(fill=batches))+theme(legend.position = "right")
	
	png(paste0(result_path,"batch_composition.png"),width = 10, height = 12, units = "in", res = 100)
	gbatch/gcelltype
	dev.off()
	# browser()
	kfile=paste0("kbet_",algorithm,"_4.rds")
	#browser()
	if(!file.exists(paste0(result_path,kfile))){
		require('FNN')
		com=run_kbet_for_each_celltype4(datasets.integrated, 
																	reduction = reduction,
																	celltype_col = celltype_col,
																	batch_col = batch_col,
																	method = method,
																	algorithm=algorithm)
		saveRDS(com, paste0(result_path,kfile))
	}else
	{
		com <- readRDS(paste0(result_path,kfile))
	}
	
	kfile=paste0("kbet_",algorithm,"_2.rds")
	#browser()
	if(!file.exists(paste0(result_path,kfile))){
	  require('FNN')
	  com2=run_kbet_for_each_celltype2(datasets.integrated, 
	                                  reduction = reduction,
	                                  celltype_col = celltype_col,
	                                  batch_col = batch_col,
	                                  method = method,
	                                  algorithm=algorithm)
	  saveRDS(com2, paste0(result_path,kfile))
	}else
	{
	  com2 <- readRDS(paste0(result_path,kfile))
	}
	
	#png(paste0(result_path,paste0("kbet_",algorithm,".png")), width=28, height=4, units="in",res=150)
	com[[2]]<-com[[2]]+ggtitle(method)
	#dev.off()
	
	#browser()
	if(!file.exists(paste0(result_path,"ari_seq_2.rds"))){
		ari_seq_celltype <- c()
		ari_seq_batch <- c()
		for(resolution in seq(0.1,1.0,by=0.1)){
			datasets.integrated <- FindClusters(datasets.integrated, resolution = resolution)
			ari_seq_celltype <- c(ari_seq_celltype,adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[celltype_col]][[1]]))
			ari_seq_batch <- c(ari_seq_batch,adjustedRandIndex(x=as.integer(datasets.integrated$seurat_clusters),y=datasets.integrated[[batch_col]][[1]]))
		}
		ari_seq_df<- data.frame(ari_celltype=ari_seq_celltype,ari_batch=ari_seq_batch)
	#   ncenters <- length(unique(celltypes))
	#   res_km<-kmeans(datasets.integrated@reductions[[reduction]]@cell.embeddings,
	#          centers = ncenters, iter.max = 30)
	#   datasets.integrated$kmclusters <- res_km$cluster
	#   ari_seq_celltype <- adjustedRandIndex(x=as.integer(datasets.integrated$kmclusters),y=datasets.integrated[[celltype_col]][[1]])
	#   ari_seq_batch <- adjustedRandIndex(x=as.integer(datasets.integrated$kmclusters),y=datasets.integrated[[batch_col]][[1]])
	# 	ari_seq_df<-data.frame(ari_celltype=ari_seq_celltype,ari_batch=ari_seq_batch)
	  saveRDS(ari_seq_df, paste0(result_path,"ari_seq_2.rds"))
	}
	else{
		ari_seq_df <- readRDS(paste0(result_path,"ari_seq_2.rds"))
	}
	if(downsample){
		ncells=dim(datasets.integrated)[2]
		datasets.integrated=datasets.integrated[,sample(1:ncells,3e4)]
	}
	# browser()
	max.k=300
	mm <- max.k-MixingMetric(datasets.integrated, grouping.var=batch_col, reduction=reduction, dims=dims, max.k=max.k)
	# sil <- silhouette(as.numeric(as.factor(datasets.integrated[[celltype_col]][[1]])),dist = dist(datasets.integrated@reductions[[reduction]]@cell.embeddings))
	sil <- silhouette(as.numeric(as.factor(datasets.integrated[[celltype_col]][[1]])),dist = dist(datasets.integrated@reductions$umap@cell.embeddings))
	# plot alignment
	library(RColorBrewer)
	# browser()
	mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(batch_type))
	png(paste0(result_path,"batch_umap.png"),width=8, height=8, units="in",res=150)
	pbatch <- DimPlot(datasets.integrated, reduction = "umap", group.by = batch_col) +
		scale_color_manual(values=mycolors)+theme_bw(base_size=30)+ 
		NoLegend() + ggtitle(method)+
		theme(axis.title = element_blank(),legend.text = element_text(size = 40,face="bold"),title = element_text(size = 70, face="bold"),
					axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5),
					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
		xlim(-limmax,limmax)+ylim(-limmax,limmax)
	plot(pbatch)
	dev.off()
	png(paste0(result_path,"celltype_umap.png"),width=8, height=8, units="in",res=150)
	pcelltype <- DimPlot(datasets.integrated, reduction = "umap", group.by = celltype_col)+
		theme_bw(base_size=30)+NoLegend()+ ggtitle(method)+
		theme(axis.title = element_blank(),legend.text = element_text(size = 40,face="bold"),title = element_text(size = 70, face="bold"),
					axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust = 0.5),
					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
		xlim(-limmax,limmax)+ylim(-limmax,limmax)
	plot(pcelltype)
	dev.off()
	
	if(run_dge)
	{
		if(!file.exists(paste0(result_path,"dge_celltype_2.rds")))
		{
			run_dge_analyse(res_path = result_path,"integrated_data_2.rds",
											batch_col = batch_col,celltype_col = celltype_col)
		}
	}
	
	return(list(median(mm),sil[,3],ari_seq_df,pbatch,pcelltype,com,com2))
}

run_dge_analyse <- function(data.integrated,
														batch1,celltype1,batch2=NULL,celltype2=NULL,
														celltype_col="celltype",batch_col="dataset",
														method="seurat")
{
	# browser()
	Idents(data.integrated) <- celltype_col
	data.integrated$celltype.cond <- paste(Idents(data.integrated), data.integrated@meta.data[[batch_col]], sep = "_")
	Idents(data.integrated) <- "celltype.cond"
	if(!is.null(batch2)){
		cc1=paste(celltype1,batch1,sep = "_")
		cc2=paste(celltype1,batch2,sep = "_")
	}
	
	if(!is.null(celltype2)){
		cc1=paste(celltype1,batch1,sep = "_")
		cc2=paste(celltype2,batch1,sep = "_")
	}
	if((!is.null(batch2)) && (!is.null(celltype2))){
		cc1=paste(celltype1,batch1,sep = "_")
		cc2=paste(celltype2,batch2,sep = "_")
	}
	nm1<-length(WhichCells(data.integrated,idents = cc1))
	nm2<-length(WhichCells(data.integrated,idents = cc2))
	print(c(nm1,nm2))
	br<-NULL
	# browser()
	if(nm1>50 && nm2>50){
		br <- FindMarkers(data.integrated, ident.1 = cc1, ident.2 = cc2, 
										 slot = "data", 
											logfc.threshold = 0.,
											test.use = "wilcox", verbose = FALSE, 
											min.cells.feature = 1,
											min.cells.group = 1,
										 min.pct = 0.1)
		# browser()
		qtest<-createData(method = "test",pval = br$p_val)
		p <- ggplot(qtest,
			aes(x=expected,y=observed,group=Method))
		p<-p+geom_ribbon(aes(ymin=clower,ymax=cupper),fill = "grey70")+
			geom_line(aes(colour=Method),size=10)+
			theme_bw(base_size=100)+theme(legend.position = "none")+ylab("Observed")+xlab("Expected")+
			theme(legend.title = element_blank())+theme(legend.text=element_text(size=70))

		print(dim(br))
		br$method=method
		br$ctype=celltype1
	}
	return(br)
}


plot_metric2 <- function(result_path,res_list){
	ari<-c()
	rr<-c()
	lrr<-c()
	gd<-list()
	methods<-names(res_list)
	i=0
	for(res in res_list){
		i=i+1
		ari <- c(ari,res[[3]])
		rr<-c(rr,res[[6]][[1]][,3])
		lrr<-c(lrr,length(res[[6]][[1]][,3]))
		print(methods[i])
		# gd[[methods[i]]]=res[[6]][[2]]
	}
	ari_data<- data.frame(cari=ari,group=rep(methods, each=10))
	kbet_data<- data.frame(mkb=rr,group=rep(methods,lrr))
	# # kBET test metric
	# png(paste0(result_path,"kbet_cvae_celltype.png"),width=28, height = 28, units="in", res=100)
	# plot_grid(plotlist = gd, labels = "auto", ncol = 1)
	# dev.off()
	png(paste0(result_path,"kbet_cvae.png"),width=4, height = 12, units="in", res=100)
	boxplot(mkb~group, data=kbet_data, ylab="kBET (acceptance rate)", las=2)
	dev.off()
	
}
plot_metric <- function(scxx_cvae_res,
												scxx_cvae2_res,
												scxx_cvae3_1_res,
												scxx_cvae3_2_res,
												scanorama_res,
												desc_res,
												seurat_res,
												mnn_res,
												result_path){
	gb<-list(scxx_cvae_res[[4]],scxx_cvae2_res[[4]],scxx_cvae3_1_res[[4]],scxx_cvae3_2_res[[4]],
					 scanorama_res[[4]],desc_res[[4]],seurat_res[[4]],mnn_res[[4]])
	gc<-list(scxx_cvae_res[[5]],scxx_cvae2_res[[5]],scxx_cvae3_1_res[[5]],scxx_cvae3_2_res[[5]],
					 scanorama_res[[5]],desc_res[[5]],seurat_res[[5]],mnn_res[[5]])
	ga<-list(scxx_cvae_res[[4]],scxx_cvae2_res[[4]],scxx_cvae3_1_res[[4]],scxx_cvae3_2_res[[4]],
					 scanorama_res[[4]],desc_res[[4]],seurat_res[[4]],mnn_res[[4]],
					 scxx_cvae_res[[5]],scxx_cvae2_res[[5]],scxx_cvae3_1_res[[5]],scxx_cvae3_2_res[[5]],
					 scanorama_res[[5]],desc_res[[5]],seurat_res[[5]],mnn_res[[5]])
	ari<-c()
	sil<-c()
	mm<-c()
	tt<-c()
	i=0
	for(res in list(scxx_cvae_res,scxx_cvae2_res,scxx_cvae3_1_res,scxx_cvae3_2_res,
									scanorama_res,desc_res,seurat_res,mnn_res)){
		i=i+1
		mm <- c(mm,res[[1]])
		sil <- c(sil,res[[2]])
		ari <- c(ari,res[[3]])
		tt<- c(tt,length(res[[2]]))
	}
	ari_data<- data.frame(cari=ari,group=rep(c("cvae","cvae2","cvae3_1","cvae3_2",
																						 "scanorama","desc","seurat","mnn"),
																					 each=10))
	sil_data<- data.frame(csil=sil,group=rep(c("cvae","cvae2","cvae3_1","cvae3_2",
																						 "scanorama","desc","seurat","mnn"),
																					 tt))
	# plot grid
	#result_path
	png(paste0(result_path,"dataset_umap.png"),width=15, height = 15, units="in", res=100)
	plot_grid(plotlist=gb, labels ="auto")
	dev.off()
	png(paste0(result_path,"celltype_umap.png"),width=15, height = 15, units="in", res=100)
	plot_grid(plotlist=gc, labels ="auto")
	dev.off()
	png(paste0(result_path,"cell_alignment_umap.png"),width=28, height = 8, units="in", res=100)
	plot_grid(plotlist = ga, labels = "auto", nrow=2)
	dev.off()
	
	# plot metric
	png(paste0(result_path,"ari.png"),width=4, height = 4, units="in", res=100)
	boxplot(cari~group,data = ari_data, ylab="Adjusted Rand Index", las=2)
	dev.off()
	# Compute or Extract Silhouette Information
	png(paste0(result_path,"silhouette.png"),width=4, height = 4, units="in", res=100)
	boxplot(csil~group,data=sil_data, ylab="Silhouette Information",las=2)
	dev.off()
	png(paste0(result_path,"mixingmetric.png"),width=4, height = 4, units="in", res=100)
	barplot(mm, ylab="Mixing Metric",names.arg=c("cvae","cvae2","cvae3_1","cvae3_2","scanorama","desc","seurat","mnn"), las=2)
	dev.off()
}

kbet_for_one_celltype<-function(rt, 
																subemb, 
																subbat,
																algorithm="cover_tree"){
	k0=floor(rt*length(subbat))
	if(k0<5){
		return(c(rt,0.05,0.05,0.05))
	}
	knn <- get.knn(subemb, k=k0, algorithm = algorithm)
	be<-kBET(df=subemb, batch = subbat, knn=knn, k0=k0, plot = FALSE)
	return(c(rt,
					 be$summary$kBET.observed[2],
					 be$summary$kBET.observed[3],
					 be$summary$kBET.observed[4]))
}
run_kbet_for_each_celltype2<- function(datasets.integrated,
																			reduction,
																			celltype_col,
																			batch_col,
																			method,
																			algorithm="cover_tree",
																			downsample=FALSE){
	# browser()
	
	if(downsample){
		ncells=dim(datasets.integrated)[2]
		datasets.integrated=datasets.integrated[,sample(1:ncells,4e3)]
	}
	embeddings <- datasets.integrated@reductions[[reduction]]@cell.embeddings
	celltypes <- as.factor(datasets.integrated@meta.data[[celltype_col]])
	batches <- datasets.integrated@meta.data[[batch_col]]
	gr<-c()
	xma<-c()
	for(f in levels(celltypes)) {
		ind= (f==celltypes)
		#ind[sample(which(ind==TRUE),round(0.9*length(which(ind==TRUE))))]=FALSE
		subemb = embeddings[ind,]
		subbat = batches[ind]
		# do foreach here.
		xm <- foreach(rt=seq(0.05, 0.25, by=0.05), 
									.combine = 'rbind',
									.export = "kbet_for_one_celltype",
									.packages = 'FNN') %dopar% kbet_for_one_celltype(rt,subemb,subbat,algorithm=algorithm)
		gr  <- c(gr,rep(f, dim(xm)[1]))
		xma <- rbind(xma,xm)
	}
	
	kb=data.frame(nb=100*xma[,1],ymax=1-xma[,2],ymedian=1-xma[,3],ymin=1-xma[,4],celltype=gr)
	h <- ggplot(kb, aes(x=nb, y=ymedian, group=celltype))
	pk <- h +
		geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70") +
		geom_line(aes(y = ymedian))+
	  theme_bw(base_size = 30)+
		facet_grid(. ~ celltype) + 
		xlab("Neighborhood size (% sample size)") + 
		ylab("KBET") + 
	  ylim(0,1)+
		labs(title = method)
	return(list(kb,pk))
}

run_kbet_for_each_celltype3<- function(datasets.integrated,
																			 reduction,
																			 celltype_col,
																			 batch_col,
																			 method,
																			 algorithm="cover_tree",
																			 downsample=FALSE){
	# browser()
	downsample=F
	if(downsample){
		ncells=dim(datasets.integrated)[2]
		datasets.integrated=datasets.integrated[,sample(1:ncells,4e3)]
	}
	embeddings <- datasets.integrated@reductions[[reduction]]@cell.embeddings
	celltypes <- as.factor(datasets.integrated@meta.data[[celltype_col]])
	batches <- datasets.integrated@meta.data[[batch_col]]
	
	xma <- foreach(rt=seq(0.05, 0.25, by=0.05), 
								.combine = 'rbind',
								.export = "kbet_for_one_celltype",
								.packages = 'FNN') %dopar% kbet_for_one_celltype(rt,embeddings,batches,algorithm=algorithm)
	
	
	kb=data.frame(nb=100*xma[,1],ymax=1-xma[,2],ymedian=1-xma[,3],ymin=1-xma[,4])
	h <- ggplot(kb, aes(x=nb, y=ymedian))
	pk <- h +
		geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70") +
		geom_line(aes(y = ymedian))+
		facet_grid(. ~ celltype) + 
		xlab("Neighborhood size (% smaple size)") + 
		ylab("KBET (acceptance rate)") + 
		labs(title = method)
	return(list(kb,pk))
}

run_kbet_for_each_celltype4<- function(datasets.integrated,
																			 reduction,
																			 celltype_col,
																			 batch_col,
																			 method,
																			 algorithm="cover_tree",
																			 downsample=FALSE){
	# browser()
	# downsample=T
	if(downsample){
		ncells=dim(datasets.integrated)[2]
		datasets.integrated=datasets.integrated[,sample(1:ncells,4e3)]
	}
	embeddings <- datasets.integrated@reductions[[reduction]]@cell.embeddings
	celltypes <- as.factor(datasets.integrated@meta.data[[celltype_col]])
	batches <- datasets.integrated@meta.data[[batch_col]]
	gr<-c()
	xma<-c()
	for(f in levels(celltypes)) {
		ind= (f==celltypes)
		subemb = embeddings[ind,]
		subbat = batches[ind]
		# do foreach here.
		xm <- foreach(rt=seq(0.05, 0.25, by=0.05), 
									.combine = 'rbind',
									.export = "kbet_for_one_celltype",
									.packages = 'FNN') %dopar% kbet_for_one_celltype(rt,subemb,subbat,algorithm=algorithm)
		gr  <- c(gr,rep(f, dim(xm)[1]))
		xma <- rbind(xma,xm)
	}
	colnames(xma)<-c("size","ymin","ymedian","ymax")
	kb<-aggregate(cbind(ymin,ymedian,ymax)~size,xma,median)
	kb=data.frame(nb=100*kb[,1],ymax=1-kb[,2],ymedian=1-kb[,3],ymin=1-kb[,4])
	h <- ggplot(kb, aes(x=nb, y=ymedian))
	pk <- h +
		geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70") +
		geom_line(aes(y = ymedian))+
		xlab("Neighborhood size (% smaple size)") + 
		ylab("KBET (acceptance rate)") + 
		labs(title = method)
	return(list(kb,pk))
}
# run_kbet_for_each_celltype<- function(datasets.integrated,
# 																			reduction,
# 																			celltype_col,
# 																			batch_col){
# 	embeddings <- datasets.integrated@reductions[[reduction]]@cell.embeddings
# 	celltypes <- as.factor(datasets.integrated@meta.data[[celltype_col]])
# 	batches <- datasets.integrated@meta.data[[batch_col]]
# 	require("FNN")
# 	rej_rate=c()
# 	for(f in levels(celltypes)){
# 		ind= (f==celltypes)
# 		subemb = embeddings[ind,]
# 		subbat = batches[ind]
# 		k0=floor(rt*length(subbat))
# 		if(k0<10)	next
# 		knn <- get.knn(subemb, k=k0, algorithm = 'cover_tree')
# 		be<-kBET(df=subemb, batch = subbat, knn=knn, k0=k0, plot = FALSE)
# 		rej_rate=c(rej_rate,be$stats$kBET.observed)
# 	}
# 	return(rej_rate)
# }

h5ad_to_seurat<-function(input_file="output.h5ad",result_path,
																	annotation_file,col.names=NA,order=NA,sep="\t"){
	# browser()
	if(!file.exists(paste0(result_path,"integrated_data.rds"))){
		datasets.integrated <- ReadH5AD(paste0(result_path,input_file),assay = "integrated")
		if(is.na(col.names)){
			metadata <- read.csv(annotation_file,
													 header = TRUE,sep = sep, stringsAsFactors = FALSE, row.names = 1)
		}else{
			metadata <- read.csv(annotation_file,
												 header = FALSE,sep = sep, stringsAsFactors = FALSE, row.names = 1,
												 col.names = col.names)
		}
		
		if(is.na(order)){
			cell_ind<-rownames(datasets.integrated@meta.data)
			if(all(rownames(datasets.integrated@meta.data)==rownames(metadata))){
				datasets.integrated@meta.data=metadata
			}else{
				datasets.integrated@meta.data=metadata[cell_ind,,drop=FALSE]
			}
		}else{
			datasets.integrated@meta.data=metadata
		}
		
		saveRDS(datasets.integrated, paste0(result_path,"integrated_data.rds"))
	}
}
scxx_to_seurat_ctr_stim<-function(){
	input_file="output.h5ad"
	result_path="./results/ctr_stim_pbmc/scxx_cvae/"
	annotation_file="./data/ctr_stim_pbmc/GSE96583/kang_pbmc_batch2_filtered/annotation.tsv"
	datasets.integrated <- ReadH5AD(paste0(result_path,input_file),assay = "integrated")
	metadata <- read.csv(annotation_file,
											 header = FALSE,sep = "\t", stringsAsFactors = FALSE, row.names = 1,
											 col.names = c("index","batch","tsne1","tsne2","ind","stim","cluster","celltype","multiplets"))
	
	if(all(rownames(datasets.integrated@meta.data)==rownames(metadata))){
		datasets.integrated@meta.data=metadata
		saveRDS(datasets.integrated, paste0(result_path,"integrated_data.rds"))
	}else{
		print("Row names don't match.")
	}
}

scxx_to_seurat_bipolar<-function(){
	input_file="output.h5ad"
	result_path="./results/bipolar/scxx_cvae/"
	infile_b1 <- "./data/bipolar/bipolarumicounts_batch1_wo_bc7/"
	infile_b2 <- "./data/bipolar/bipolarumicounts_batch2_wo_bc2/"
	metadata.b1 <- read.csv(paste0(infile_b1,"annotation.tsv"),header = FALSE,sep = "\t", row.names = 1, col.names = c("barcode","celltype","cluster","replicate_id","batch_name"),stringsAsFactors = FALSE)
	metadata.b2 <- read.csv(paste0(infile_b2,"annotation.tsv"),header = FALSE,sep = "\t", row.names = 1, col.names = c("barcode","celltype","cluster","replicate_id","batch_name"),stringsAsFactors = FALSE)
	datasets.integrated <- ReadH5AD(paste0(result_path,input_file),assay = "integrated")
	metadata <- rbind(metadata.b1, metadata.b2)
	if(all(rownames(datasets.integrated@meta.data)==rownames(metadata))){
		datasets.integrated@meta.data=metadata
		saveRDS(datasets.integrated, paste0(result_path,"integrated_data.rds"))
	}else{
		print("Row names don't match.")
	}
}


source("./code/base/qqplot.R")
plot_volcano <- function(br.all, method, ctype){
	br.use <-br.all[br.all$method==method & br.all$ctype==ctype,]
	pvol<-EnhancedVolcano(br.use,
												lab = rownames(br.use),
												x = 'avg_logFC',
												y = 'p_val_adj',
												xlim = c(-4, 4),
												pCutoff = 1e-20,
												FCcutoff = 1.5,
												title = paste0(method," : ", substr(ctype,1,6)),
												xlab = bquote(~Log[2]~'fold change'),
												ylab = bquote(~-Log[10]~adjusted~italic(P)),
												legend = c("NS","Log2 FC","p_val_adj","p_val_adj & Log2 FC")
												)
	zero_ind<-br.use$p_val<1e-300
	if(any(zero_ind)){
		br.use[zero_ind,]$p_val=1e-300
	}
	#pptl <- qqunif.plot(br.use$p_val)
	df <- createData(method, br.use$p_val)
	df$celltype=ctype
	return(list(pvol,df))
}

## partII: trajectory
run_3steps <- function(cds_obj, result_path,  title, method,start_cells="1h", num_dim = 10){
	cds_obj <- preprocess_cds(cds_obj, num_dim = num_dim, norm_method = "none")
	cds_obj <- cluster_cells(cds_obj, reduction_method="UMAP")
	cds_obj <- learn_graph(cds_obj)
	root_pr_nodes = get_earliest_principal_node(cds_obj,startcell = start_cells )
	cds_obj <- order_cells(cds_obj, root_pr_nodes = root_pr_nodes)
	p1<-plot_cells(cds_obj, label_groups_by_cluster=FALSE, 
								 color_cells_by = "batch", label_cell_groups = FALSE, 
								 trajectory_graph_segment_size=1.5,
								 cell_size = 6, group_label_size = 10)+ggtitle(method)+
								theme(legend.direction = "horizontal")+
								theme(legend.title = element_blank(),
											legend.text = element_text(size = 60),
											legend.position = c(0.4,0.5))
	# pseudotime on reduced dimension
	reducedDim(cds_obj,"UMAP") <- reducedDim(cds_obj,"RD")
	cds_obj <- cluster_cells(cds_obj, reduction_method="UMAP")
	cds_obj <- learn_graph(cds_obj)
	root_pr_nodes = get_earliest_principal_node(cds_obj,startcell = start_cells )
	cds_obj <- order_cells(cds_obj, root_pr_nodes = root_pr_nodes)
	saveRDS(cds_obj, file = paste0(result_path,title,"_cds_pt_on_dr.rds"))
	return(p1)
}
run_monocle3<- function(cds.integrated,result_path,t1,t2,method,
												split_by="condition",num_dim=20,result_file="traj_expdata.png"){
	# browser()
	cds_1 <- cds.integrated[,cds.integrated@colData[[split_by]]==t1]
	cds_2 <- cds.integrated[,!(cds.integrated@colData[[split_by]]==t1)]
	pl_1<-run_3steps(cds_1,result_path=result_path, start_cells="1h", num_dim=num_dim, title=t1,method=method)
	pl_2<-run_3steps(cds_2,result_path=result_path, start_cells="1h", num_dim=num_dim, title=t2,method=method)
	png(paste0(result_path,result_file),width=12, height = 12, units="in", res=100)
	g<-plot_grid(plotlist = list(pl_1,pl_2), nrow = 2)
	print(g)
	dev.off()
	return(list(pl_1,pl_2))
}
run_monocle3_integrated<- function(cds.integrated,result_path,t1,t2,method,start_cells="4W",title="LPS & PAM",
												split_by="condition",num_dim=20,result_file="traj_expdata_integrated.png"){
	#browser()
	g<-run_3steps(cds.integrated,result_path=result_path, 
								start_cells=start_cells, num_dim=num_dim, title=title,method=method)
	png(paste0(result_path,result_file),width=12, height = 12, units="in", res=100)
	print(g)
	dev.off()
	return(g)
}
# plot trajectory for visualization.
plot_trajectory <- function(datasets.integrated, result_path, t1="PAM",t2="LPS",method="", col.name=NA, reduction="scxx", pseudotime_on_dr=FALSE){
	# trajectory using scale data
  # browser()
	datasets.integrated <- run_umap_on_seurat(datasets.integrated, result_path, reduction=reduction)
	expression_data = datasets.integrated@assays$integrated@scale.data
	cell_metadata = datasets.integrated@meta.data
	gene_metadata = datasets.integrated@assays$integrated@meta.features
	if(is.na(col.name)){
		gene_metadata["gene_short_name"]=rownames(gene_metadata)
	}
	else{
		gene_metadata["gene_short_name"]=gene_metadata[col.name]
	}
	cds.integrated <- new_cell_data_set(expression_data = expression_data,
																			cell_metadata = cell_metadata,
																			gene_metadata = gene_metadata)
	reducedDim(cds.integrated,"RD")=datasets.integrated@reductions[[reduction]]@cell.embeddings
	reducedDim(cds.integrated,"UMAP")=datasets.integrated@reductions[["umap"]]@cell.embeddings
	saveRDS(cds.integrated, file = paste0(result_path,"cds_expressdata.rds"))
	g<-run_monocle3(cds.integrated, result_path=result_path,t1=t1,t2=t2,method=method,result_file = "traj_corrected_data.png")
	return(g)
}
plot_trajectory2 <- function(datasets.integrated, result_path, t1="PAM",t2="LPS",method="", col.name=NA, reduction="scxx", pseudotime_on_dr=FALSE, assay="integrated"){
	# browser()
	#datasets.integrated <- run_umap_on_seurat(datasets.integrated, result_path, reduction=reduction)
	expression_data = datasets.integrated[[assay]]@scale.data
	cell_metadata = datasets.integrated@meta.data
	gene_metadata = datasets.integrated[[assay]]@meta.features
	if(is.na(col.name)){
		gene_metadata["gene_short_name"]=rownames(gene_metadata)
	}
	else{
		gene_metadata["gene_short_name"]=gene_metadata[col.name]
	}
	cds.integrated <- new_cell_data_set(expression_data = expression_data,
																			cell_metadata = cell_metadata,
																			gene_metadata = gene_metadata)
	reducedDim(cds.integrated,"RD")=datasets.integrated@reductions[[reduction]]@cell.embeddings
	reducedDim(cds.integrated,"UMAP")=datasets.integrated@reductions[["umap"]]@cell.embeddings
	saveRDS(cds.integrated, file = paste0(result_path,"cds_expressdata.rds"))
	g<-run_monocle3_integrated(cds.integrated, result_path=result_path,t1=t1,t2=t2,method=method,result_file = "traj_corrected_data.png")
	return(g)
}


plot_trajectory_germline <- function(datasets.integrated, result_path, t1="F",t2="M",
																		 start_cells="1h", title="Li & Guo",
																		 method="", col.name=NA, reduction="scxx", pseudotime_on_dr=FALSE){
	# browser()
	expression_data = datasets.integrated@assays$integrated@scale.data
	cell_metadata = datasets.integrated@meta.data
	gene_metadata = datasets.integrated@assays$integrated@meta.features
	if(is.na(col.name)){
		gene_metadata["gene_short_name"]=rownames(gene_metadata)
	}
	else{
		gene_metadata["gene_short_name"]=gene_metadata[col.name]
	}
	cds.integrated <- new_cell_data_set(expression_data = expression_data,
																			cell_metadata = cell_metadata,
																			gene_metadata = gene_metadata)
	reducedDim(cds.integrated,"RD")=datasets.integrated@reductions[[reduction]]@cell.embeddings
	reducedDim(cds.integrated,"UMAP")=datasets.integrated@reductions[["umap"]]@cell.embeddings
	saveRDS(cds.integrated, file = paste0(result_path,"cds_expressdata.rds"))
	g<-run_monocle3_integrated(cds.integrated, result_path=result_path,t1=t1,t2=t2,method=method,result_file = "traj_corrected_data.png")
	return(g)
}

run_umap_on_seurat <- function(datasets.integrated, reduction="scxx",min.dist=0.1){
	dims=1:16
	datasets.integrated <- RunUMAP(datasets.integrated, reduction = reduction, dims = dims, min.dist = min.dist, n_threads=10)
	datasets.integrated <- RunTSNE(datasets.integrated, reduction = reduction, dims=dims, check_duplicates = FALSE)
	return(datasets.integrated)
}

## the one with the largest frequency in the start cell stage treated as the starter
get_earliest_principal_node <- function(input,startcell){
	cell_ids       <- which(colData(input)[, "week"] == startcell)
	closest_vertex <- input@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
	closest_vertex <- as.matrix(closest_vertex[colnames(input), ])
	root_pr_nodes  <- igraph::V(principal_graph(input)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
	return(root_pr_nodes)
}

writeSeurat2mtx<-function(obj,data_path){
	if(!dir.exists(data_path)){
		dir.create(data_path,recursive = TRUE)
	}
	mat <- obj@assays$RNA@counts
	mat@x <- floor(mat@x)
	meta.data <- obj@meta.data
	meta.features <- obj@assays$RNA@meta.features
	write.table(meta.data,file = paste0(data_path,"barcodes.tsv"),quote = FALSE,sep = "\t",col.names=FALSE)
	write.table(meta.features,file=paste0(data_path,"genes.tsv"),quote = FALSE,sep = "",col.names=FALSE)
	writeMM(mat,file=paste0(data_path,"matrix.mtx"))
}

celltypePrediction<-function(adata,ref,query){
	Idents(adata)<-"X_batch"
	adata_ref<-adata[,WhichCells(adata,idents = ref)]
	adata_query<-adata[,WhichCells(adata,idents = query)]
	Idents(adata_ref)<-"IdentChar"
	M<-c()
	for(ct in levels(adata$IdentChar)){
		#print(ct)
		cell_ind<-WhichCells(adata_ref,idents = ct)
		Reductions(adata_ref,slot = "scxx")
		centroid<-apply(adata_ref@reductions[["scxx"]]@cell.embeddings[cell_ind,],2,median)
		M<-cbind(M,centroid)
	}
	colnames(M)<-levels(adata$IdentChar)
	Q<-adata_query@reductions[["scxx"]]@cell.embeddings
	D<-matrix(0,dim(Q)[1],dim(M)[2])
	for(i in 1:dim(Q)[1]){
		for(j in 1:dim(M)[2]){
			D[i,j]<-dist(rbind(Q[i,],M[,j]), method = "euclidean")
		}
	}
	rownames(D)<-rownames(Q)
	colnames(D)<-colnames(M)
	res<-list()
	res$prediction<-colnames(D)[apply(D, 1, which.min)]
	res$score<-exp(-(apply(D, 1, min))**2/2)
	res$celltype<- as.character(adata_query$IdentChar)
	res<-as.data.frame(res)
	return(res)
}

rescale<-function(x){
	return((x-min(x))/(max(x)-min(x)))
}

