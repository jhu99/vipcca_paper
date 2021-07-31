rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)

### create gene activity matrix
peaks<-Read10X_h5("./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
# activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "./data/atacseq/10x_pbmc/Homo_sapiens.GRCh37.82.gtf",
# 																						seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)
# saveRDS(activity.matrix,"./data/atacseq/10x_pbmc/geneActMat.rds")
activity.matrix<-readRDS("./data/atacseq/10x_pbmc/geneActMat.rds")
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_singlecell.csv", sep = ",", header = TRUE, row.names = 1,
									 stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
pbmc.atac$tech <- "atac"
#
# ## data preprocessing
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
pbmc.rna<-readRDS("./data/atacseq/10x_pbmc/pbmc_10k_v3.rds")
pbmc.rna$tech<-"rna"

## integrated atac-seq and rna-seq with cca
genelist<-intersect(rownames(pbmc.atac),VariableFeatures(pbmc.rna))
genelist<-genelist[rowSums(pbmc.atac@assays$ACTIVITY@counts)[genelist]>0]
VariableFeatures(pbmc.rna) <- genelist
VariableFeatures(pbmc.atac) <- genelist
min(rowSums(pbmc.atac@assays$ACTIVITY@counts)[genelist])

reference.list <- list(batch1=pbmc.rna, batch2=pbmc.atac)
b.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:16)
b.integrated <- IntegrateData(anchorset = b.anchors, dims = 1:16)
b.integrated <- ScaleData(b.integrated, verbose = FALSE)
b.integrated <- RunPCA(b.integrated, npcs = 16, verbose = FALSE)
b.integrated <- RunUMAP(b.integrated, dims = 1:16)
saveRDS(b.integrated,"./results/atacseq/10x_pbmc/b.integrated.rds")
b.integrated<-readRDS("./results/atacseq/10x_pbmc/b.integrated.rds")
b.integrated$tech <- as.factor(b.integrated$tech)
levels(b.integrated$tech)<-c("ATAC","RNA")
b.integrated$celltype<- as.factor(b.integrated$celltype)
levels(b.integrated$celltype)<-c("B cell progenitor","CD14+ Monocytes", "CD16+ Monocytes", "CD4 Memory","CD4 Naive", "CD8 effector","CD8 Naive","Dendritic cell", "Double negative T cell","NK cell",  "pDC", "Platelets", "Pre-B cell")

p1<-DimPlot(b.integrated, group.by = "tech")+ggtitle("Seurat V3")+labs(tag="a")+
	theme(legend.position = c(0.7,0.9),axis.title = element_blank(),
				plot.tag = element_text(size = 100),text  = element_text(size=60),
				legend.text = element_text(size = 70,face = "bold"),title = element_text(size=70),
				axis.text = element_blank(),axis.ticks = element_blank())+
	guides(colour = guide_legend(override.aes = list(size=15)))

p2<-DimPlot(b.integrated, group.by = "celltype",split.by = "tech")+
	theme(axis.title = element_blank(),text  = element_text(size=60),
				legend.text = element_text(size = 70),
				axis.text = element_blank(),axis.ticks = element_blank())+
	guides(colour = guide_legend(override.aes = list(size=15)))+ggtitle("")
png("./results/atacseq/10x_pbmc/row1.png",width = 50,height = 15, res = 100, units = "in")
plot_grid(p1,p2,nrow = 1,rel_widths = c(2,3))
dev.off()

## kbet
source("./code/run/evaluation_function.R")
b.integrated@meta.data[["is_selected"]]="all_cells_types"
kfile="kbet_seurat_4.rds"
celltype_col="is_selected"
batch_col="tech"
reduction="pca"
dims=1:16
method="seurat"
result_path="./results/atacseq/10x_pbmc/Seurat/"
if(!file.exists(paste0(result_path,kfile,sep=""))){
	require('FNN')
	#browser()
	com_s=run_kbet_for_each_celltype4(b.integrated,
																	reduction = reduction,
																	celltype_col = celltype_col,
																	batch_col = batch_col,
																	method = method,
																	downsample = TRUE)
	saveRDS(com_s, paste0(result_path,kfile,sep=""))
}else{
	com_s <- readRDS(paste0(result_path,kfile,sep=""))
}
## mixing metric
max.k=300
mm_s <- max.k-MixingMetric(b.integrated, grouping.var=batch_col, reduction=reduction, dims=dims, max.k=max.k)
median(mm_s)



# integrated with vipcca
library(SeuratDisk)
fd<-"./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16"
Convert(paste(fd,"output.h5ad",sep = "/"), dest = "h5seurat", overwrite = TRUE)
adata<-LoadH5Seurat(paste(fd,"output.h5seurat",sep = "/"))
VariableFeatures(adata)<-rownames(adata)
#adata<-ReadH5AD(paste(fd,"output.h5ad",sep = "/"))
ann<-read.csv(paste(fd,"annotation.csv",sep = "/"),row.names = 1)
adata<-RenameCells(adata,new.names=rownames(ann))
ann_atac<-read.csv(paste(fd,"ann_atac.csv",sep = "/"),row.names = 1,stringsAsFactors = F)
ann_atac[ann_atac$pred_max_score<=0,'celltype']=NA
ann[rownames(ann_atac),'celltype']=ann_atac$celltype
adata<-AddMetaData(adata,metadata = ann)
adata$tech <- as.factor(adata$tech)
levels(adata$tech)<-c("ATAC","RNA")
adata<-RunUMAP(adata, reduction = "scxx", dims = 1:16)


p1<-DimPlot(b.integrated, group.by = "tech")+ggtitle("Seurat V3")+
	theme(axis.title = element_blank(),legend.position = c(0.2,0.5),
				plot.tag = element_text(size = 40),text  = element_text(size=28),
				legend.text = element_text(size = 25,face="bold"),title = element_text(size=28),
				axis.text = element_blank(),axis.ticks = element_blank())+
	guides(colour = guide_legend(override.aes = list(size=10)))
p2<-DimPlot(b.integrated, group.by = "celltype",split.by = "tech")+
	theme(axis.title = element_blank(),text  = element_text(size=28),
				legend.text = element_text(size = 25,face = "bold"),legend.direction = "horizontal",
				axis.text = element_blank(),axis.ticks = element_blank())+
	guides(colour = guide_legend(override.aes = list(size=10),nrow=5))+ggtitle("")
legend1<-cowplot::get_legend(p1)
legend2<-cowplot::get_legend(p2)
p1<-p1+NoLegend()
p2<-p2+NoLegend()
row1<-plot_grid(p1,p2,nrow = 1,rel_widths = c(1,2))+labs(tag = "a")+theme(plot.tag = element_text(size=40),plot.margin = unit(c(1,0,0,1),"cm"))
p3<-DimPlot(adata, group.by = "tech")+ggtitle("VIPCCA")+
	theme(axis.title = element_blank(),
				plot.tag = element_text(size = 40),text  = element_text(size=28),
				legend.text = element_text(size = 28,face = "bold"),title = element_text(size=28),
				axis.text = element_blank(),axis.ticks = element_blank())+NoLegend()
## do it latter
p4<-DimPlot(adata, group.by = "celltype",split.by = "tech")+
	theme(axis.title = element_blank(),legend.title = element_blank(),legend.direction = "horizontal",
				plot.tag = element_text(size = 40),text  = element_text(size=28),
				legend.text = element_text(size = 28,face = "bold"),title = element_text(size=28),
				axis.text = element_blank(),axis.ticks = element_blank())+NoLegend()
	
row2<-plot_grid(p3,p4,nrow = 1,rel_widths = c(1,2))+labs(tag = "b")+theme(plot.tag = element_text(size=40),plot.margin = unit(c(1,0,0,1),"cm"))
row3<-plot_grid(legend1,legend2,rel_widths = c(1,2))
#part1<-{(p1+p2+plot_layout(widths = c(1,2)))/(p3+p4+plot_layout(widths = c(1,2)))}+plot_layout(heights = c(2,3))
part1<-plot_grid(plotlist = list(row1,row2,row3),ncol = 1,rel_heights = c(3,3,1))
png("./results/atacseq/10x_pbmc/part1.png", width = 20, height = 14,units = "in",res = 100)
plot(part1)
dev.off()

## kbet on vipcca
source("./code/run/evaluation_function.R")
adata@meta.data[["is_selected"]]="all_cells_types"
kfile="kbet_seurat_4.rds"
celltype_col="is_selected"
batch_col="tech"
reduction="scxx"
dims=1:16
method="vipcca"
result_path="./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16/"
if(!file.exists(paste0(result_path,kfile,sep=""))){
	require('FNN')
	#browser()
	com_v=run_kbet_for_each_celltype4(adata,
																	reduction = reduction,
																	celltype_col = celltype_col,
																	batch_col = batch_col,
																	method = method,
																	downsample = TRUE)
	saveRDS(com_v, paste0(result_path,kfile,sep=""))
}else{
	com_v <- readRDS(paste0(result_path,kfile,sep=""))
}
## mixing metric
max.k=300
mm_v <- max.k-MixingMetric(adata, grouping.var=batch_col, reduction=reduction, dims=dims, max.k=max.k)
median(mm_v)
# mixing metric on vipcca
colors=c("cyan","chocolate","coral","chartreuse",
				 "burlywood","blue","orchid","deeppink1")

## kbet curves
df_com<-rbind(com_v[[1]],com_s[[1]])
df_com$method <- as.factor(rep(c("VIPCCA","Seurat V3"),each=5))
png("./results/atacseq/10x_pbmc/kbet.png",width = 15,height = 8,units = "in",res = 100)
gkbet<-ggplot(df_com, aes(x=nb,y=ymedian))+
	geom_line(aes(group=method,color=method),size=5)+
	scale_color_manual(values=colors[c(7,8)])+theme(text = element_text(size=40,face = "bold"),legend.position = c(0.2,0.5),
													panel.grid.major = element_blank(),legend.key.size = unit(2,"cm"),
													panel.grid.minor = element_blank(),legend.title = element_blank(),
													panel.background = element_blank(),
													axis.line = element_line(colour = "black"))+
													xlab("Neighborhood size (%)")+ylab("kBET")
gkbet
dev.off()


library(SingleCellExperiment)
library(scmap)
data.rna<-SingleCellExperiment(assays=list(logcounts=as.matrix(GetAssayData(pbmc.rna))),
															 colData=pbmc.rna@meta.data)
rowData(data.rna)$feature_symbol <- rownames(data.rna)
rowData(data.rna)$scmap_features<-FALSE
rowData(data.rna)[genelist,"scmap_features"]=TRUE
data.rna<-indexCluster(data.rna,cluster_col = "celltype")
heatmap(as.matrix(metadata(data.rna)$scmap_cluster_index))
data.activity<-SingleCellExperiment(assays=list(logcounts=as.matrix(GetAssayData(pbmc.atac,assay = "ACTIVITY"))))
rowData(data.activity)$feature_symbol <- rownames(data.activity)
scmapCluster_results <- scmapCluster(
	projection = data.activity, 
	index_list = list(
		yan = metadata(data.rna)$scmap_cluster_index
	)
)
table(scmapCluster_results$scmap_cluster_labs)
####scmap-cell
set.seed(1)
data.rna<-indexCell(data.rna)
str(metadata(data.rna)$scmap_cell_index$subclusters)
dim(metadata(data.rna)$scmap_cell_index$subcentroids[[1]])

scmapCell_results <- scmapCell(
	data.activity, 
	list(
		yan = metadata(data.rna)$scmap_cell_index
	)
)
scmapCell_results$yan$similarities[,1:3]
scmapCell_clusters <- scmapCell2Cluster(
	scmapCell_results, 
	list(
		as.character(colData(data.rna)$cell_type)
	)
)
table(scmapCell_clusters$combined_labs)


###########

celltype.predictions<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
celltype.predictions[celltype.predictions$prediction.score.max<0.5,1]="Unknown"
celltype.predictions<-celltype.predictions[,c(1,15)]
ann_atac[is.na(ann_atac$celltype),"celltype"]="Unknown"
table(ann_atac$celltype)
table(celltype.predictions$predicted.id)
table(ann_atac$celltype,celltype.predictions$predicted.id)

xname<-rownames(ann_atac)
yname<-rownames(celltype.predictions)
df_jaccard <-data.frame()
for(x in levels(factor(ann_atac$celltype))){
	for(y in levels(factor(celltype.predictions$predicted.id))){
		xsubset<- xname[ann_atac$celltype==x]
		ysubset<- yname[celltype.predictions$predicted.id==y]
		xlist<-list(vipcca=x,seurat=y,
								Jaccard=(1.0*length(intersect(xsubset,ysubset)))/length(union(xsubset,ysubset)),
								Overlap=length(intersect(xsubset,ysubset)))
		df_jaccard <- rbind(df_jaccard,as.data.frame(xlist))
	}
}
df_jaccard<- df_jaccard[df_jaccard$Overlap>0,]
df_jaccard[which(df_jaccard$vipcca=="pre-B cell"),1]="Pre-B cell"
df_jaccard[which(df_jaccard$seurat=="pre-B cell"),2]="Pre-B cell"


p5<-ggplot(df_jaccard, aes(x=vipcca,y=seurat,color=Jaccard,size=Overlap))+
	geom_point()+labs(tag = "c")+ scale_size(range = c(0, 30))+
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5),axis.text = element_text(size = 28,face = "bold"),
				plot.tag = element_text(size = 40,face="bold"),text = element_text(size = 28),
				axis.title = element_text(size = 28,face = "bold"),panel.grid.major = element_blank(),
				plot.margin = unit(c(1,0,0,1),"cm"),
				panel.grid.minor=element_blank(),panel.background = element_blank())+
	scale_color_gradient(low="blue", high="red")+labs(x="VIPCCA",y="Seurat V3")

png("./results/atacseq/10x_pbmc/part2.png", width = 20, height = 9,units = "in",res = 100)
plot(p5)
dev.off()

#
adata.atac <- subset(adata, subset=tech=="ATAC")
sum(rownames(celltype.predictions)==colnames(adata.atac))
adata.atac@meta.data$predicted.id<-celltype.predictions$predicted.id
adata.atac@meta.data$is_unknown_by_vipcca=is.na(adata.atac$celltype)
adata.atac@meta.data$is_unknown_by_seurat=(adata.atac$predicted.id=="Unknown")
x<-pbmc.atac@meta.data[,c("duplicate","chimeric","unmapped","TSS_fragments","blacklist_region_fragments","mitochondrial")]
adata.atac<-AddMetaData(adata.atac,metadata = x)
Idents(adata.atac)<-"is_unknown_by_vipcca"
p6<-FeaturePlot(adata.atac,reduction = "umap", pt.size = 3,
								features  = c("is_unknown_by_seurat","is_unknown_by_vipcca","duplicate","chimeric","unmapped","TSS_fragments"),combine = T) &
	theme(axis.title = element_blank(),text = element_text(size=24),axis.text = element_blank(),axis.ticks = element_blank())
p6[[1]]<-p6[[1]]+ggtitle("Unknown (Seurat V3)")
p6[[2]]<-p6[[2]]+ggtitle("Unknown (VIPCCA)")
p6[[3]]<-p6[[3]]+ggtitle("Duplicate")
p6[[4]]<-p6[[4]]+ggtitle("Chimeric")
p6[[5]]<-p6[[5]]+ggtitle("Unmapped")
p6[[6]]<-p6[[6]]+ggtitle("TSS fragments")
part3<-	plot_grid(p6,label_x = 0, label_y = 1,hjust = 0.5, vjust = 0.5,nrow = 1)+labs(tag="d")+theme(plot.tag.position = "topleft",plot.margin = unit(c(1,0,0,1),"cm"),
																plot.tag = element_text(size = 40,face = "bold"))
png("./results/atacseq/10x_pbmc/part3.png", width = 12, height = 7,units = "in",res = 100)
plot(part3)
dev.off()

df_unknown_seurat<-adata.atac@meta.data[adata.atac@meta.data[,"is_unknown_by_seurat"],
										 c("duplicate","chimeric","unmapped","TSS_fragments")]
df_unknown_vipcca<-adata.atac@meta.data[adata.atac@meta.data[,"is_unknown_by_vipcca"],
																				c("duplicate","chimeric","unmapped","TSS_fragments")]
df_unknown_counts<-data.frame(
meanV=c(colMeans(df_unknown_seurat),colMeans(df_unknown_vipcca)),
methods=rep(c("Seurat V3","VIPCCA"),each=4),
items=rep(c("Duplicate","Chimeric","Unmapped","TSS fragments"),2)
)
df_unknown_counts_1<-df_unknown_counts[c(1,4,5,8),]
df_unknown_counts_2<-df_unknown_counts[c(2,3,6,7),]



pp1<-ggplot(df_unknown_counts_1,aes(x=items, y=meanV,fill=methods))+
	geom_bar(stat = "identity",position=position_dodge())+
	scale_fill_manual(values = colors[c(7,8)])+
	labs(tag="e")+
	theme(text = element_text(size=35,face = "bold"),axis.title.x = element_blank(),axis.title.y = element_blank(),
				plot.tag=element_text(size=40),legend.position = c(0.8,0.9),legend.title = element_blank(),
				panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
				panel.background=element_blank(),plot.margin = unit(c(1,0,0,1),"cm"))
pp2<-ggplot(df_unknown_counts_2,aes(x=items, y=meanV,fill=methods))+
	geom_bar(stat = "identity",position=position_dodge())+
	scale_fill_manual(values = colors[c(7,8)])+
	labs(tag="f")+
	theme(text = element_text(size=35,face = "bold"),axis.title = element_blank(),
				panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
				panel.background=element_blank(),plot.margin = unit(c(1,0,0,1),"cm"),
				plot.tag=element_text(size=40))+NoLegend()
part4<-plot_grid(pp1,pp2,ncol=1)
row3<-plot_grid(part3,part4,nrow = 1)
png("./results/atacseq/10x_pbmc/part4.png", width = 8, height = 9,units = "in",res = 100)
pp1/pp2
dev.off()

png("./results/atacseq/10x_pbmc/all-in-one.png",width = 22,height = 32,units = "in",res = 100)
plot_grid(plotlist = list(part1,p5,row3),ncol = 1,rel_heights = c(14,9,9))
#{p5/part3/(pp1+pp2)}+plot_layout(ncol = 1,heights = c(10,10,10))
dev.off()





png("./results/atacseq/10x_pbmc/fig5.png",width = 30,height = 50,units = "in",res = 100)
#plot_grid(p1,p2,rel_widths=c(1,2))
{(p1+p2+plot_layout(widths = c(1,2)))/(p3+p4+plot_layout(widths = c(1,2)))/(p5)/p6/(pp1+pp2)}+plot_layout(heights = c(1,1,2,3,2))
dev.off()

# write celltype to files
seurat.labels<-adata.atac@meta.data[,13,drop=FALSE]
seurat.labels<-seurat.labels[seurat.labels$predicted.id!="Unknown",,drop=F]
seurat.labels$predicted.id<-gsub(' ','-',seurat.labels$predicted.id)
vipcca.labels<-adata.atac@meta.data[,3,drop=FALSE]
vipcca.labels[is.na(vipcca.labels$celltype),]="Unknown"
vipcca.labels<-vipcca.labels[vipcca.labels$celltype!="Unknown",,drop=F]
vipcca.labels$celltype<-gsub(' ','-',vipcca.labels$celltype)
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
t1<-t1[,4:15]
df1<-expand.grid(marker=rownames(t1),celltype=colnames(t1))
df1$coverage <- as.vector(as.matrix(t1))
df1$method <- "Seurat V3"

t2<-read.csv(f2,sep = "\t",quote = "\'")
rownames(t2)<-markers2
t2<-t2[markers1,]
t2<-t2[,4:15]
df2<-expand.grid(marker=rownames(t2),celltype=colnames(t2))
df2$coverage <- as.vector(as.matrix(t2))
df2$method <- "VIPCCA"
df_coverage<-rbind(df1,df2)
df_coverage<-df_coverage[with(df_coverage,order(marker,celltype)),]

library(ggplot2)
png("./results/atacseq/10x_pbmc/coverage.png",width = 1000,height = 1000)
p7<-ggplot(df_coverage, aes(x=celltype,y=coverage,fill=method))+
	geom_bar(stat = "identity",position=position_dodge())+
	facet_grid(marker~.,scales = "free_y")+theme_bw(base_size = 30)+
	ylab("Average read coverage of chromatin accessibility (RPKM)")+theme(legend.position = "top")+theme(legend.title = element_blank())
p7
dev.off()

# genome track
### atac track plots
#library(Signac)
# Error: Dependency package(s) 'Rsamtools','ggseqlogo','ggbio',
# 'biovizBase','AnnotationFilter','fastmatch','lsa','RcppRoll','qvalue' not available
# fastmatch','lsa','RcppRoll','qvalue' 
rm(list = ls())
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


seurat.labels<-adata.atac@meta.data[,13,drop=FALSE]
vipcca.labels<-adata.atac@meta.data[,3,drop=FALSE]
vipcca.labels[is.na(vipcca.labels$celltype),]="Unknown"

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
a1<-a[1:2]
a2<-a[5:8]
#
subtitle=c("GNLY","MS4A1","BCL11B","LYZ","CD8A","CD3F","CD4","HLA-DRA")
devtools::load_all("../signac/")
undebug(CoveragePlot)
p8<-CoveragePlot(
	object = pbmc,
	region = a[1:2],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "celltype",
	annotation = FALSE,
	peaks = FALSE,
	tag="a",
	title="VIPCCA"
)&theme(text = element_text(size = 15,face = "bold"))
p9<-CoveragePlot(
	object = pbmc,
	region = a[1:2],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	annotation=TRUE,
	peaks=TRUE,
	links=TRUE,
	assay = "peaks",
	group.by = "predicted.id",
	tag="b",
	title="Seurat V3"
)&theme(text = element_text(size = 15,face = "bold"))

p10<-wrap_plots(p8,p9,ncol = 1,heights = c(4,6))
png("./results/atacseq/10x_pbmc/genome_track_1.png",width = 10, height = 10, units = "in", res = 100)
plot(p10)
dev.off()


p11<-CoveragePlot(
	object = pbmc,
	region = a[3:4],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "celltype",
	annotation = FALSE,
	peaks = FALSE,
	tag="a",
	title="VIPCCA"
)&theme(text = element_text(size = 15,face = "bold"))
p12<-CoveragePlot(
	object = pbmc,
	region = a[3:4],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	annotation=TRUE,
	peaks=TRUE,
	links=TRUE,
	assay = "peaks",
	group.by = "predicted.id",
	tag="b",
	title="Seurat V3"
)&theme(text = element_text(size = 15,face = "bold"))

png("./results/atacseq/10x_pbmc/genome_track_2.png",width = 12, height = 12, units = "in", res = 100)
wrap_plots(p11,p12,ncol = 1,heights = c(4,6))
dev.off()


p13<-CoveragePlot(
	object = pbmc,
	region = a[5:6],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "celltype",
	annotation = FALSE,
	peaks = FALSE,
	tag="a",
	title="VIPCCA"
)&theme(text = element_text(size = 15,face = "bold"))
p14<-CoveragePlot(
	object = pbmc,
	region = a[5:6],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	annotation=TRUE,
	peaks=TRUE,
	links=TRUE,
	assay = "peaks",
	group.by = "predicted.id",
	tag="b",
	title="Seurat V3"
)&theme(text = element_text(size = 15,face = "bold"))

png("./results/atacseq/10x_pbmc/genome_track_3.png",width = 12, height = 12, units = "in", res = 100)
wrap_plots(p13,p14,ncol = 1,heights = c(4,6))
dev.off()


p15<-CoveragePlot(
	object = pbmc,
	region = a[7:8],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	assay = "peaks",
	group.by = "celltype",
	annotation = FALSE,
	peaks = FALSE,
	tag="a",
	title="VIPCCA"
)&theme(text = element_text(size = 15,face = "bold"))
p16<-CoveragePlot(
	object = pbmc,
	region = a[7:8],
	sep = c(":", "-"),
	extend.upstream = 10000,
	extend.downstream = 2000,
	nrow = 1,
	annotation=TRUE,
	peaks=TRUE,
	links=TRUE,
	assay = "peaks",
	group.by = "predicted.id",
	tag="b",
	title="Seurat V3"
)&theme(text = element_text(size = 15,face = "bold"))

png("./results/atacseq/10x_pbmc/genome_track_4.png",width = 10, height = 10, units = "in", res = 100)
wrap_plots(p15,p16,ncol = 1,heights = c(4,6))
dev.off()


###
#download bulk ATAC count file
url<-"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74912&format=file&file=GSE74912%5FATACseq%5FAll%5FCounts%2Etxt%2Egz"
download.file(url,destfile = "./data/atacseq/10x_pbmc/GSE74912_ATACseq_All_Counts.txt.gz")
counts<-read.delim("./data/atacseq/10x_pbmc/GSE74912_ATACseq_All_Counts.txt")
# download bulk ATAC BAM file
# url for ucsc genome hub track
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=https://s3-us-west-1.amazonaws.com/chang-public-data/2016_NatGen_ATAC-AML/hub.txt
	
dir.create("./data/atacseq/bulk_atac_pbmc/",showWarnings = F)
download.file("https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/CD4-9.merge.s20.w150sw.bw",
							destfile = "./data/atacseq/bulk_atac_pbmc/CD4.bw")
download.file("https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/CD8-10.merge.s20.w150sw.bw",
							destfile = "./data/atacseq/bulk_atac_pbmc/CD8.bw")
download.file("https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/Nkcell-11.merge.s20.w150sw.bw",
							destfile = "./data/atacseq/bulk_atac_pbmc/NKcell.bw")
download.file("https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/Bcell-13.merge.s20.w150sw.bw",
							destfile = "./data/atacseq/bulk_atac_pbmc/Bcell.bw")
download.file("https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/Mono-7.merge.s20.w150sw.bw",
							destfile = "./data/atacseq/bulk_atac_pbmc/mono.bw")

# 
t1 <- read.table("./vipcca_bam/scores_per_transcript_genome.tab",sep = "\t",na.strings = "nan",quote = "'",comment.char = "",header = T)
t2 <- read.table("./seurat_bam/scores_per_transcript_genome.tab",sep = "\t",na.strings = "nan",quote = "'",comment.char = "",header = T)
t3 <- read.table("./data/atacseq/bulk_atac_pbmc/scores_per_transcript_genome.tab",sep = "\t",na.strings = "nan",quote = "'",comment.char = "",header = T)

t1[,2]<- as.character(t1[,2])
t1[,3]<- as.character(t1[,3])
xname1<-t1[,1:3]
xname1<-apply(xname1,MARGIN = 1,function(x) paste0(x,collapse =":"))
t1<-t1[,-1:-3]
rownames(t1)<-xname1

t2[,2]<- as.character(t2[,2])
t2[,3]<- as.character(t2[,3])
xname2<-t2[,1:3]
xname2<-apply(xname2,MARGIN = 1,function(x) paste0(x,collapse =":"))
t2<-t2[,-1:-3]
rownames(t2)<-xname2

t3[,2]<- as.character(t3[,2])
t3[,3]<- as.character(t3[,3])
xname3<-t3[,1:3]
xname3<-apply(xname3,MARGIN = 1,function(x) paste0(x,collapse =":"))
t3<-t3[,-1:-3]
rownames(t3)<-xname3
commonname<-intersect(intersect(xname1,xname2),xname3)
colnames(t1)
t1<-t1[commonname,]
t2<-t2[commonname,]
t3<-t3[commonname,]
df_cor <- data.frame()


cxs<-colnames(t1)[c(1,4,2,11,6,9,7,5,8)]

for(cx in cxs){
	for(cy in colnames(t3)){
		print(c(cx,cy))
		
		df_cor<-rbind(df_cor,
									as.data.frame(list(corvipcca=unname(cor.test(t1[,cx],t3[,cy],method="pearson")$estimate),
											 corseurat=unname(cor.test(t2[,cx],t3[,cy],method = "pearson")$estimate),
											 singlecell=cx,
											 bulk=cy)))
	}
}
df_cor$sign = 0
df_cor[c(1,6,12,17,23,28,34,40,45),5]=1
sum(df_cor$corvipcca*df_cor$sign)/9
sum(df_cor$corseurat*df_cor$sign)/9

df_cor$singlecell
library(tidyverse)

p11<-ggplot(df_cor, aes(x=singlecell,y=bulk, fill=corvipcca))+
	geom_tile()+ggtitle("VIPCCA (pearson: 0.80)")+labs(x="scATAC-seq",y="Bulk chromatin accessibility")+
	scale_fill_gradient(limits=c(0.0,1.0),name="correlation",low = "white",high = "red")+
	theme(text = element_text(size = 40,face = "bold"),axis.title = element_text(size = 40))+NoLegend()
p11<-p11+geom_rect(color="black",size=1,mapping = aes(xmin=c(0.5),xmax=c(2.5),	ymin=c(0.5),ymax=c(1.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(2.5),xmax=c(4.5),	ymin=c(1.5),ymax=c(2.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(4.5),xmax=c(6.5),	ymin=c(2.5),ymax=c(3.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(6.5),xmax=c(7.5),	ymin=c(3.5),ymax=c(4.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(7.5),xmax=c(9.5),	ymin=c(4.5),ymax=c(5.5)))

p12<-ggplot(df_cor, aes(x=singlecell,y=bulk, fill=corseurat))+
	geom_tile()+ggtitle("Seurat V3 (pearson : 0.80)")+labs(x="scATAC-seq",y="bulk chromatin accessibility")+
	scale_fill_gradient(limits=c(0.0,1.0),name="correlation",low = "white",high = "red")+
	theme(text = element_text(size = 40,face="bold"),axis.title = element_text(size = 40),
				legend.direction = "horizontal",legend.key.size = unit(3,"cm") )
p12<-p12+geom_rect(color="black",size=1,mapping = aes(xmin=c(0.5),xmax=c(2.5),	ymin=c(0.5),ymax=c(1.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(2.5),xmax=c(4.5),	ymin=c(1.5),ymax=c(2.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(4.5),xmax=c(6.5),	ymin=c(2.5),ymax=c(3.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(6.5),xmax=c(7.5),	ymin=c(3.5),ymax=c(4.5)))+
	geom_rect(color="black",size=1,mapping = aes(xmin=c(7.5),xmax=c(9.5),	ymin=c(4.5),ymax=c(5.5)))
legend_correlation<-cowplot::get_legend(p12)
p12<-p12+NoLegend()
gpearson<-plot_grid(p11,p12,nrow = 1,rel_heights = c(1,1))
png("./results/atacseq/10x_pbmc/pearson_correlation.png",width = 15,height = 5,units = "in",res = 100)
gpearson
dev.off()
png("./results/atacseq/10x_pbmc/sp.png",width = 30,height = 30,units = "in",res = 100)
plot_grid(plotlist=list(gkbet,gpearson,legend_correlation),ncol = 1,labels = c("a","b"),rel_heights = c(4,4,1),
					label_size = 70,label_x = 0,label_y = 1,vjust = 0.5,hjust = 0.5)+
	theme(plot.margin = margin(1,0,0,1,unit = "cm"))
dev.off()

# coembedding ... lsi
celltype.predictions<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
celltype.predictions[celltype.predictions$prediction.score.max<0.5,1]="unknown"
celltype.predictions<-celltype.predictions[,c(1,15)]

