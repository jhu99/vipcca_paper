rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)

# ### create gene activity matrix
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
# # 
# # ## data preprocessing
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

pbmc.atac<-readRDS("./data/atacseq/10x_pbmc/pbmc_atac.rds")
# RNA datasets
# download.file("https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1","./data/atacseq/10x_pbmc/pbmc_10k_v3.rds")
pbmc.rna<-readRDS("./data/atacseq/10x_pbmc/pbmc_10k_v3.rds")
pbmc.rna$tech<-"rna"
# DimPlot(pbmc.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
# DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")

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

transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
																				reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype,
																		 weight.reduction = pbmc.atac[["lsi"]])

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
pbmc.atac$predicted.id[pbmc.atac$prediction.score.max<=0.5]<-"unknown"
p1 <- DimPlot(pbmc.atac, reduction = "umap",group.by = "predicted.id",label = TRUE, repel = TRUE) + ggtitle("scATAC-seq")+NoLegend()
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
p1+p2
# hist(pbmc.atac$prediction.score.max)
# abline(v = 0.5, col = "red")
# table(pbmc.atac$prediction.score.max > 0.5)
# saveRDS(celltype.predictions,"./results/atacseq/10x_pbmc/celltype.prediction.rds")
# 
# View the prediction
pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.5)
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = levels(pbmc.rna))  # to make the colors match
p3 <- DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") +
	NoLegend() + scale_colour_hue(drop = FALSE)
p4 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") +
	NoLegend()
remove(pbmc.atac.filtered)
CombinePlots(plots = list(p1,p2,p3,p4))
# 
# # Co-embedding
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]])
pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc.rna, y = pbmc.atac)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
p5 <- DimPlot(coembed, group.by = "tech")+ggtitle("Seurat V3")+theme_bw(base_size = 30)+theme(legend.position = c(0.6,0.8))
p6 <- DimPlot(coembed, group.by = "celltype",split.by="tech",combine = T)+theme_bw(base_size = 30)
#+ggtitle("Seurat V3")+theme_bw(base_size = 30)+NoLegend()
p5+p6


# pbmc.rna$seurat_clusters
# f<-read.table("./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_cDNA/features.tsv.gz")
#rm(list = ls())
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
# celltype.predictions<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
# celltype.predictions[celltype.predictions$prediction.score.max<0.5,1]="unknown"
# celltype.predictions<-celltype.predictions[,c(1,15)]
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

adata<-ReadH5AD(paste(fd,"output.h5ad",sep = "/"))
ann<-read.csv(paste(fd,"annotation.csv",sep = "/"),row.names = 1)
ann_atac<-read.csv(paste(fd,"ann_atac.csv",sep = "/"),row.names = 1,stringsAsFactors = F)
ann_atac[ann_atac$pred_max_score<1e-300,'celltype']='unknown'
ann_atac["CTTCTAATCACTACCC-1",'celltype']="unknown"
ann[rownames(ann_atac),'celltype']=ann_atac$celltype

#
adata<-AddMetaData(adata,metadata = ann)
adata<-RunUMAP(adata, reduction = "scxx", dims = 1:16)
p7<-DimPlot(adata, group.by = "tech")+ggtitle("VIPCCA")+theme_bw(base_size = 30)+theme(legend.position = c(0.6,0.8))
p8<-DimPlot(adata, group.by = "celltype",split.by = "tech")+ggtitle("VIPCCA")+theme_bw(base_size = 30)
legend2<-cowplot::get_legend(p8)


#datasets<-SplitObject(adata,split.by = "tech")
adata.rna<-subset(adata, subset = tech=="rna")
adata.atac<-subset(adata, subset = tech=="atac")

adata.atac@meta.data$predicted.id<-pbmc.atac$predicted.id

for(ct in levels(adata.atac$celltype)){
	print(ct)
}

pbmc.atac$predicted.id<-as.factor(pbmc.atac$predicted.id)
levels(pbmc.atac$predicted.id)<-c(levels(pbmc.atac$predicted.id),"Platelets")
df_heatmap<-as.data.frame(table(adata.atac$celltype,pbmc.atac$predicted.id))
df_heatmap$Freq<-log1p(df_heatmap$Freq)
df_heatmap<-df_heatmap[with(df_heatmap,order(Var1,Var2)),]
levels(df_heatmap$Var1)<-c("Bpro","CD14M","CD16","CD4M","CD4N","CD8E",
													 "CD8N","DC","DNT","NK","pDC","PLA","Bpre","unknown")
levels(df_heatmap$Var2)<-c("Bpro","CD14M","CD16","CD4M","CD4N","CD8E",
													 "CD8N","DC","DNT","NK","pDC","Bpre","unknown","PLA")
df_heatmap$Var1<-as.character(df_heatmap$Var1)
df_heatmap$Var2<-as.character(df_heatmap$Var2)
df_heatmap<-df_heatmap[with(df_heatmap,order(Var1,Var2)),]
p9<-ggplot(df_heatmap,aes(x = Var2,y=Var1,fill=Freq))+geom_tile()+
	xlab("Seurat V3")+ylab("VIPCCA")+
	labs(fill="Scaled Freq")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p9
adata.atac@meta.data$is_unknown_by_vipcca=(adata.atac$celltype=="unknown")
adata.atac@meta.data$is_unknown_by_seurat=(adata.atac$predicted.id=="unknown")

x<-pbmc.atac@meta.data[,c("duplicate","chimeric","unmapped","TSS_fragments","blacklist_region_fragments","mitochondrial")]
adata.atac<-AddMetaData(adata.atac,metadata = x)

p10<-FeaturePlot(adata.atac,reduction = "umap", features  = c("is_unknown_by_seurat","is_unknown_by_vipcca","duplicate","chimeric","unmapped","TSS_fragments"),combine = T)
png("./results/atacseq/10x_pbmc/fig_sx1.png",width = 1000,height = 500)
p9+p10
dev.off()
DefaultAssay(pbmc.atac)<-"ACTIVITY"
pbmc.atac.coembed=pbmc.atac
pbmc.atac.coembed@reductions$umap<-NULL
pbmc.atac.coembed[["umap"]]<-coembed[,colnames(pbmc.atac.coembed)]@reductions[["umap"]]
p11<-FeaturePlot(pbmc.atac.coembed,reduction = "umap", features = c("GNLY","MS4A1","BCL11B","LYZ","CD8A","CD3E","CD4","HLA-DRA"),combine = F)
pbmc.atac.coembed[["umap"]]<-adata.atac[,colnames(pbmc.atac.coembed)]@reductions[["umap"]]
p12<-FeaturePlot(pbmc.atac.coembed,reduction = "umap", features = c("GNLY","MS4A1","BCL11B","LYZ","CD8A","CD3E","CD4","HLA-DRA"),combine = F)
title_text<-c("GNLY (NK)","MS4A1 (Bpro,Bpre)","BCL11B (CD4N,CD8E,CD8N,CD4M)","LYZ (CD14M)",
				 "CD8A (CD8E, CD8N)","CD3E (CD4N, CD8E, NK, CD8N)","CD4 (CD4N, CD4M, CD14M, CD16)","HLA-DRA (CD16, CD14M)")
for(i in 1:8){
	p11[[i]]<-p11[[i]]+ggtitle(title_text[i])
	p12[[i]]<-p12[[i]]+ggtitle(title_text[i])
}

png("./results/atacseq/10x_pbmc/fig_sx2.png",width=2000,height = 2000)
plot_grid(plotlist = c(p11,p12),ncol = 4)
dev.off()


seurat.labels<-adata.atac@meta.data[,15,drop=FALSE]
seurat.labels<-seurat.labels[seurat.labels$predicted.id!="unknown",,drop=F]
seurat.labels$predicted.id<-sub(' ','-',seurat.labels$predicted.id)

vipcca.labels<-adata.atac@meta.data[,3,drop=FALSE]
vipcca.labels<-vipcca.labels[vipcca.labels$celltype!="unknown",,drop=F]
vipcca.labels$celltype<-sub(' ','-',vipcca.labels$celltype)
write.table(vipcca.labels,
						"./data/atacseq/10x_pbmc/vipcca.barcode.celltype.txt",
						quote = FALSE, col.names = FALSE,sep = "\t")
write.table(seurat.labels,
						"./data/atacseq/10x_pbmc/seurat.barcode.celltype.txt",
						quote = FALSE,col.names = FALSE,sep="\t")

# accessibility enrichment around marker genes
markers1<-c("GNLY","MS4A1","BCL11B","LYZ","CD8A","CD3E","CD4","HLA-DRA")
markers2<-c("FCGR3A","ITGB1","IRF7","MS4A1","CD3E","CD4","LYZ","BCL11B","GNLY","CD8A","HLA-DRA")
f1 <- "./seurat_bam/scores_per_transcript.tab"
f2 <- "./vipcca_bam/scores_per_transcript.tab"
t1<-read.csv(f1, sep = "\t",quote = "\'")
rownames(t1)<-markers2
t1<-t1[markers1,]
t1$PLA=0
t1<-t1[,4:16]
df1<-expand.grid(marker=rownames(t1),celltype=colnames(t1))
df1$coverage <- as.vector(as.matrix(t1))
df1$method <- "Seurat V3"

t2<-read.csv(f2,sep = "\t",quote = "\'")
rownames(t2)<-markers2
t2<-t2[markers1,]
t2$pDC=0
t2$DC=0
t2$PLA=0
t2<-t2[,4:16]
df2<-expand.grid(marker=rownames(t2),celltype=colnames(t2))
df2$coverage <- as.vector(as.matrix(t2))
df2$method <- "VIPCCA"

df_coverage<-rbind(df1,df2)
df_coverage<-df_coverage[with(df_coverage,order(marker,celltype)),]


library(ggplot2)
png("./results/atacseq/10x_pbmc/coverage.png",width = 2000,height = 1500)
p13<-ggplot(df_coverage, aes(x=celltype,y=coverage,fill=method))+
	geom_bar(stat = "identity",position=position_dodge())+
	facet_grid(marker~.,scales = "free_y")+theme_bw(base_size = 30)+
	ylab("Average read coverage of chromatin accessibility (RPKM)")+theme(legend.position = "top")+theme(legend.title = element_blank())
p13
dev.off()


pr1<-plot_grid(p5,p6,p7,p8,ncol = 2,
							 rel_widths = c(1,2),labels = c("a","b","c","d"),label_size = 30)
png("./results/atacseq/10x_pbmc/all_plot.png",width = 2000, height = 2000)
plot_grid(pr1,p13,ncol = 1,rel_heights = c(1,1),labels = c("","e"),label_size = 30)
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


seurat.labels<-adata.atac@meta.data[,15,drop=FALSE]
vipcca.labels<-adata.atac@meta.data[,3,drop=FALSE]

# seurat.labels<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
# seurat.labels[seurat.labels$prediction.score.max<0.5,1]="unknown"
# seurat.labels<-seurat.labels[,c(1,15)]
# vipcca.labels<-read.csv("./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16/annotation.csv",row.names = 1)
# vipcca.labels<-vipcca.labels[vipcca.labels$tech=="atac",c( "size_factor","celltype")]

pbmc<-subset(pbmc,cells=rownames(seurat.labels))
pbmc <- AddMetaData(object = pbmc, metadata = seurat.labels)
pbmc <- AddMetaData(object = pbmc, metadata = vipcca.labels)

#Idents(pbmc)<-"predicted.id"

a<- c("chr2:85912298-85925978","chr11:60223225-60238234","chr14:99635624-99737862",
			"chr12:69742121-69748015","chr2:87011729-87035520","chr11:118175260-118186891",
			"chr12:6896024-6929975","chr1:161501549-161602918","chr6:32397619-32414824",
			"chr11:602553-618000","chr10:33179247-33296721")
#a<-a[1:2]
a<-a[c(-8,-10,-11)]
devtools::load_all("../signac/")
p14<-CoveragePlot(
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
	tag="f",
	title="VIPCCA"
)
p15<-CoveragePlot(
	object = pbmc,
	region = a,
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "predicted.id",
	tag="g",
	title="Seurat V3"
)

png("./results/atacseq/10x_pbmc/all_plot_2.png",width = 20, height = 18, units = "in", res = 100)
p16<-wrap_plots(p14,p15,ncol = 1,heights = c(4,5))
p16
dev.off()

