rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)

#### create gene activity matrix
# peaks<-Read10X_h5("./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
# # activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "./data/atacseq/10x_pbmc/Homo_sapiens.GRCh37.82.gtf", 
# # 																						seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)
# # saveRDS(activity.matrix,"./data/atacseq/10x_pbmc/geneActMat.rds")
# activity.matrix<-readRDS("./data/atacseq/10x_pbmc/geneActMat.rds")
# pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
# pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
# meta <- read.table("./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
# 									 stringsAsFactors = FALSE)
# meta <- meta[colnames(pbmc.atac), ]
# pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
# pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
# pbmc.atac$tech <- "atac"
# 
# 
# ## data preprocessing
# DefaultAssay(pbmc.atac) <- "ACTIVITY"
# pbmc.atac <- FindVariableFeatures(pbmc.atac)
# pbmc.atac <- NormalizeData(pbmc.atac)
# pbmc.atac <- ScaleData(pbmc.atac)
# # dir.create("./data/atacseq/10x_pbmc/geneActMat/",showWarnings = FALSE)
# # writeMM(pbmc.atac@assays$ACTIVITY@counts,"./data/atacseq/10x_pbmc/geneActMat/matrix.mtx")
# # write.table(colnames(pbmc.atac),"./data/atacseq/10x_pbmc/geneActMat/barcodes.tsv",
# # 						row.names = FALSE,col.names = FALSE,quote = FALSE)
# # write.table(cbind(rownames(pbmc.atac),rownames(pbmc.atac)),"./data/atacseq/10x_pbmc/geneActMat/genes.tsv",sep = "\t",
# # 						row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# DefaultAssay(pbmc.atac) <- "ATAC"
# VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))
# pbmc.atac <- RunLSI(pbmc.atac, n = 50, scale.max = NULL)
# pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:50)
# #saveRDS(pbmc.atac,"./data/atacseq/10x_pbmc/pbmc_atac.rds")

# pbmc.atac<-readRDS("./data/atacseq/10x_pbmc/pbmc_atac.rds")
# # RNA datasets
# # download.file("https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1","./data/atacseq/10x_pbmc/pbmc_10k_v3.rds")
# pbmc.rna<-readRDS("./data/atacseq/10x_pbmc/pbmc_10k_v3.rds")
# pbmc.rna$tech<-"rna"
# # p1 <- DimPlot(pbmc.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
# # p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
# # CombinePlots(plots = list(p1, p2))
# # write to 10x_mtx
# library(Matrix)
# dir.create("./data/atacseq/10x_pbmc/ExpressMat/",showWarnings = FALSE)
# writeMM(pbmc.rna@assays$RNA@counts,"./data/atacseq/10x_pbmc/ExpressMat/matrix.mtx")
# write.table(colnames(pbmc.rna),"./data/atacseq/10x_pbmc/ExpressMat/barcodes.tsv",
# 						row.names = FALSE,col.names = FALSE,quote = FALSE)
# write.table(cbind(rownames(pbmc.rna),rownames(pbmc.rna)),"./data/atacseq/10x_pbmc/ExpressMat/genes.tsv",sep = "\t",
# 						row.names = FALSE, col.names = FALSE, quote = FALSE)
# geneset1<-rownames(pbmc.atac)
# geneset2<-VariableFeatures(pbmc.rna)
# commongenes<- intersect(geneset1,geneset2)
# write.table(commongenes,file = "./data/atacseq/10x_pbmc/genelist.txt",
# 						row.names = FALSE, col.names = FALSE, quote = FALSE)

# transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
# 																				reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
# celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype, 
# 																		 weight.reduction = pbmc.atac[["lsi"]])
# 
# pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
# 
# hist(pbmc.atac$prediction.score.max)
# abline(v = 0.5, col = "red")
# table(pbmc.atac$prediction.score.max > 0.5)
# saveRDS(celltype.predictions,"./results/atacseq/10x_pbmc/celltype.prediction.rds")
# 
# # View the prediction
# pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.5)
# pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = levels(pbmc.rna))  # to make the colors match
# p1 <- DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
# 	NoLegend() + scale_colour_hue(drop = FALSE)
# p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
# 	NoLegend()
# CombinePlots(plots = list(p1, p2))
# 
# # Co-embedding
# genes.use <- VariableFeatures(pbmc.rna)
# refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]
# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]])
# pbmc.atac[["RNA"]] <- imputation
# coembed <- merge(x = pbmc.rna, y = pbmc.atac)
# coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
# coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
# coembed <- RunUMAP(coembed, dims = 1:30)
# coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
# p1 <- DimPlot(coembed, group.by = "tech")
# p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
# CombinePlots(list(p1, p2))

# pbmc.rna$seurat_clusters
# f<-read.table("./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_cDNA/features.tsv.gz")
#rm(list = ls())
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
celltype.predictions<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
celltype.predictions[celltype.predictions$prediction.score.max<0.5,1]="unknown"
celltype.predictions<-celltype.predictions[,c(1,15)]
library(rlist)
fd<-"./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16"
# 
# pbmc.rna<-readRDS("./data/atacseq/10x_pbmc/pbmc_10k_v3.rds")
# pbmc.rna$tech="rna"
# ann.rna<-pbmc.rna@meta.data[,c("celltype","tech")]
# pbmc.atac<-readRDS("./data/atacseq/10x_pbmc/pbmc_atac.rds")
# ann.atac<-pbmc.atac@meta.data
# ann.atac$celltype<-"unknown"
# ann.atac<-ann.atac[,c("celltype","tech")]
# ann<-rbind(ann.rna,ann.atac)
# saveRDS(ann,"./data/atacseq/10x_pbmc/annotation.rds")
# ann<-readRDS("./data/atacseq/10x_pbmc/annotation.rds")
ga<- list()
discrepance<-c()
for(fd in dir("./results/atacseq/10x_pbmc/scxx_cvae",pattern = "cvae_5_128_16",full.names = T)){
	#if(file.exists())
	adata<-ReadH5AD(paste(fd,"output.h5ad",sep = "/"))
	ann<-read.csv(paste(fd,"annotation.csv",sep = "/"),row.names = 1)
	adata<-AddMetaData(adata,metadata = ann)
	adata<-RunUMAP(adata, reduction = "scxx", dims = 1:16)
	
	#datasets<-SplitObject(adata,split.by = "tech")
	adata.rna<-subset(adata, subset = tech=="rna")
	adata.atac<-subset(adata, subset = tech=="atac")
	adata.atac<-AddMetaData(adata.atac,metadata = celltype.predictions)
	
	for(ct in levels(adata.atac$celltype)){
		print(ct)
	}
	df_heatmap<-as.data.frame(table(adata.atac$celltype,adata.atac$predicted.id))
	df_heatmap$Freq<-log1p(df_heatmap$Freq)
	df_heatmap<-df_heatmap[with(df_heatmap,order(Var1,Var2)),]
	levels(df_heatmap$Var1)<-c("Bpro","CD14M","CD16","CD4M","CD4N","CD8E",
														 "CD8N","DC","DNT","NK","pDC","PLA","Bpre","unknown")
	levels(df_heatmap$Var2)<-c("Bpro","CD14M","CD16","CD4M","CD4N","CD8E",
														 "CD8N","DC","DNT","NK","pDC","Bpre","unknown")

	p5l1<-ggplot(df_heatmap,aes(x = Var2,y=Var1,fill=Freq))+geom_tile()+
		xlab("Seurat V3")+ylab("VIPCCA")+
		labs(fill="Scaled Freq")+
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	
	discrepance<-c(discrepance,sum(adata.atac$celltype!=adata.atac$predicted.id))
	indicator<-(adata.atac$celltype!=adata.atac$predicted.id)
	df1<-as.data.frame(table(adata.atac[,indicator]$celltype))
	df1$method="VIPCCA"
	df2<-as.data.frame(table(adata.atac[,indicator]$predicted.id))
	df2$method="Seurat V3"
	df<-rbind(df1,df2)
	#saveRDS(adata.atac,paste(fd,"atac.rds",sep = "/"))
	#adata.atac<-readRDS(paste(fd,"atac.rds",sep = "/"))
	adata.atac@meta.data$is_unknown_vipcca=(adata.atac$celltype=="unknown")
	adata.atac@meta.data$is_unknown_seurat=(adata.atac$predicted.id=="unknown")
	pbmc.atac<-readRDS("./data/atacseq/10x_pbmc/pbmc_atac.rds")
	pbmc.atac@reductions$umap@cell.embeddings<-adata.atac@reductions$umap@cell.embeddings
	FeaturePlot(pbmc.atac,reduction = "umap", features = "TSS_fragments")
	p4l1<-FeaturePlot(pbmc.atac,reduction = "umap", features  = c("duplicate","chimeric","unmapped","TSS_fragments"),combine = FALSE)
	p4l2<-FeaturePlot(adata.atac,reduction = "umap",  features  = c("is_unknown_vipcca","is_unknown_seurat"),combine = FALSE)
	DefaultAssay(pbmc.atac)<-"ACTIVITY"
	p4l3<-FeaturePlot(pbmc.atac,reduction = "umap", features = c("LYZ","GNLY","MS4A1","BCL11B","CD8A","CD3E","CD4","FCGR3A","HLA-DRA","ITGB1"),combine = FALSE)
	p4<-c(p4l2,p4l1,p4l3)
	png("./results/atacseq/10x_pbmc/feature_plots.png",width = 20,height = 20,units = "in",res = 200)
	plot_grid(plotlist = p4,ncol = 4,labels = "auto")
	dev.off()
	png("./results/atacseq/10x_pbmc/feature_plot_marker.png",width = 12,height = 8,units = "in",res = 200)
	p4l3
	dev.off()
	p5<-ggplot(data = df,aes(x=method,y=Freq,fill=Var1))+
			geom_bar(position="dodge", stat="identity")+xlab("")+
		ylab("# inconsistent prediction")+theme(legend.title =element_blank())
	legend1<-cowplot::get_legend(p5)
	p5<-p5+NoLegend()
	p6<-DimPlot(adata, reduction = "umap", group.by = "tech")+ggtitle("ATAC + RNA")
	p7<-DimPlot(adata.rna, reduction = "umap", group.by = "celltype",label = TRUE, repel = TRUE)+ggtitle("RNA")+NoLegend()
	p8<-DimPlot(adata.atac, reduction = "umap", group.by = "celltype",label = TRUE, repel = TRUE)+ggtitle("ATAC (VIPCCA)")+NoLegend()
	p9<-DimPlot(adata.atac, reduction = "umap",group.by = "predicted.id",label = TRUE, repel = TRUE)+ggtitle("ATAC (Seurat V3)")+NoLegend()
	
	# pa<-c(list(p6,p7,p8,p9,p5l1),p4)
	# #pa<-list.append(pa,legend1)
	# png(paste(fd,"paired_plot_a.png",sep="/"),width = 12, height = 8, units = "in", res = 300)
	# g<-plot_grid(plotlist = pa,nrow = 4,
	# 					labels = "auto")
	# print(g)
	# dev.off()
}
#saveRDS(discrepance,"./results/atacseq/10x_pbmc/scxx_cvae/discrepance.rds")


# accessibility enrichment around marker genes
markers<-c("LYZ","GNLY","MS4A1","BCL11B","CD8A","CD3E","CD4","FCGR3A","HLA-DRA","IRF7","ITGB1")
f1 <- "./seurat_bam/scores_per_transcript.tab"
f2 <- "./vipcca_bam/scores_per_transcript.tab"
t1<-read.csv(f1, sep = "\t",quote = "\'")
rownames(t1)<-markers
t1<-t(t1[,4:16])
ct1<-paste(rownames(t1),"S",sep = ".")
t2<-read.csv(f2,sep = "\t",quote = "\'")
rownames(t2)<-markers
t2<-t(t2[,4:17])
ct2<-paste(rownames(t2),"V",sep = ".")
library(gtools)
t<-smartbind(t1,t2,fill = 0)
rownames(t)<-c(ct1,ct2)
t<-t[order(rownames(t)),]
t<-t[1:25,]
cellgroups=rownames(t)
p10.list<-list()
df.acc.markers<-expand.grid(cellgroups,markers)
df.acc.markers$acc<-as.vector(unlist(t))

p10<-ggplot(data = df.acc.markers,aes(x=Var1,y=acc))+
	geom_bar(stat="identity")+
	facet_grid(rows  = vars(Var2),scales = "free")+
	ylab("Average normalized Accessibility")+xlab("Predicted cell types")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pr1<-plot_grid(plotlist = list(p6,p7,p8,p9),nrow=1,labels = "auto")
pr2<-plot_grid(p5l1,p10,nrow = 1,rel_widths = c(1,1),labels = c('e',"f"))
png("./results/atacseq/10x_pbmc/all_plot.png",width = 16, height = 12, units = "in", res = 300)
plot_grid(pr1,pr2,ncol = 1,rel_heights = c(1,2))
dev.off()

### atac track plots
#library(Signac)
devtools::load_all("../signac/")
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(zoo)
library(dplyr)
library(gggenes)
library(rlist)
library(cowplot)
set.seed(1234)
counts_path<-"./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
meta_path<-"./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_singlecell.csv"
fragment.path <- "./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz"
#source('./code/atacseq_scripts/Signac_function.R')

counts <- Read10X_h5(filename = counts_path)
# 
metadata <- read.csv(
	file = meta_path,
	header = TRUE,
	row.names = 1
)
chrom_assay <- CreateChromatinAssay(
	counts = counts,
	sep = c(":", "-"),
	min.cells = 10,
	genome="hg38",
	fragments = fragment.path
)
#
pbmc <- CreateSeuratObject(
	counts = chrom_assay,
	assay = "peaks",
	meta.data = metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

# gene.ranges <- genes(EnsDb.Hsapiens.v75)
# seqlevelsStyle(gene.ranges) <- 'UCSC'
# gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
# gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

#pbmc <- NucleosomeSignal(object = pbmc)
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
# pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
# 
# pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
# FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

seurat.labels<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
seurat.labels[seurat.labels$prediction.score.max<0.5,1]="unknown"
seurat.labels<-seurat.labels[,c(1,15)]
vipcca.labels<-read.csv("./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16/annotation.csv",row.names = 1)
vipcca.labels<-vipcca.labels[vipcca.labels$tech=="atac",c( "size_factor","celltype")]
pbmc<-subset(pbmc,cells=rownames(seurat.labels))
pbmc <- AddMetaData(object = pbmc, metadata = seurat.labels)
pbmc <- AddMetaData(object = pbmc, metadata = vipcca.labels)

#Idents(pbmc)<-"predicted.id"

a<- c("chr12:69742121-69748015","chr2:85912298-85925978","chr11:60223225-60238234",
			"chr14:99635624-99737862","chr2:87011729-87035520","chr11:118175260-118186891",
			"chr12:6896024-6929975","chr1:161501549-161602918","chr6:32397619-32414824",
			"chr11:602553-618000","chr10:33179247-33296721")
#a<-a[1:2]
devtools::load_all("../signac/")
p1<-CoveragePlot(
	object = pbmc,
	region = a,
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "celltype",
	annotation = FALSE,
	peaks = FALSE,
	tag="g",
	title="VIPCCA"
)
p2<-CoveragePlot(
	object = pbmc,
	region = a,
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "predicted.id",
	tag="h",
	title="Seurat V3"
)

png("./results/atacseq/10x_pbmc/all_plot_2.png",width = 16, height = 12, units = "in", res = 100)
wrap_plots(p1,p2,ncol = 1,heights = c(8,10))
dev.off()
