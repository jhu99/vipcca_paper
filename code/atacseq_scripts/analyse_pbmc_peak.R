rm(list = ls())
library(Signac)
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
url="http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz"
#download.file(url = url,destfile = "./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz")
url="http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
#download.file(url = url,destfile = "./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz.tbi")
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
#pbmc <- RunTFIDF(pbmc)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

gene.ranges <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

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
	peaks = FALSE
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
)

#p3<-CombineTracks(plotlist = list(p1,p2),heights = c(8,10))
png("./results/atacseq/10x_pbmc/accessibility_track_1.png",width = 8*12,height = 8*8,units = "in",res = 100)
p3<-wrap_plots(p1,p2,nrow=2,heights = c(8,10))
plot(p3)
dev.off()

# "./data/atacseq/10x_pbmc/atac_v1_pbmc_10k_possorted_bam.bam"
seurat.labels<-readRDS("./results/atacseq/10x_pbmc/celltype.prediction.rds")
seurat.labels[seurat.labels$prediction.score.max<0.5,1]="unknown"
seurat.labels<-seurat.labels[,1,drop=FALSE]
vipcca.labels<-read.csv("./results/atacseq/10x_pbmc/scxx_cvae/scmb_cvae_5_128_16/annotation.csv",row.names = 1)
vipcca.labels<-vipcca.labels[vipcca.labels$tech=="atac","celltype",drop=FALSE]

write.table(vipcca.labels,
			"./data/atacseq/10x_pbmc/vipcca.barcode.celltype.txt",
			quote = FALSE, col.names = FALSE,sep = "\t")
write.table(seurat.labels,
						"./data/atacseq/10x_pbmc/seurat.barcode.celltype.txt",
						quote = FALSE,col.names = FALSE,sep="\t")

# 
# bulk.atac.file<-"./data/atacseq/10x_pbmc/GSE74912_ATACseq_All_Counts.txt.gz"
# bulk.atac.count<-read.csv(bulk.atac.file,sep = "\t",stringsAsFactors = FALSE)
# paste(bulk.atac.count[1:10,1:3],sep = "-")
# 

# 0	IL7R, CCR7	Naive CD4+ T
# 1	IL7R, S100A4	Memory CD4+
# 2	CD14, LYZ	CD14+ Mono
# 3	MS4A1	B
# 4	CD8A	CD8+ T
# 5	FCGR3A, MS4A7	FCGR3A+ Mono
# 6	GNLY, NKG7	NK
# 7	FCER1A, CST3	DC
# 8	PPBP	Platelet

# cd16 monocytes markers:
# CD14 and CD16, CCR2, CD36, HLA-DR, and CD11c 
# DC markers:



# remove cd14 unreasonable signal 
# 
# IRF7: pDC
# CD29, CD31, CD36, CD41: Platelet
# HLA-DR+ cd56-: DC Human dendritic cell subsets
# CD16: CD16+ Monocytes
# CD14: CD14+ Monocytes



chrs<-c("chr12","chr2","chr11","chr14","chr2","chr11","chr12","chr1","chr6","chr11","chr10")
markers<-c("LYZ","GNLY","MS4A1","BCL11B","CD8A","CD3E","CD4","FCGR3A","HLA-DRA","IRF7","ITGB1")
start_pos<-c()
end_pos<-c()
for(mk in markers){
	marker_index <- which(gene.ranges@elementMetadata@listData$symbol== mk)
	start<-gene.ranges@ranges@start[marker_index]
	end<-start+gene.ranges@ranges@width[marker_index]
	start_pos<-c(start_pos,start-10000)
	end_pos<-c(end_pos,end+2000)
}
genebed<-data.frame(chrom=chrs,start=start_pos,end=end_pos)
write.table(genebed,file = "./data/atacseq/10x_pbmc/markers.bed",
					row.names = FALSE,sep = "\t",col.names = FALSE,
					quote = FALSE)
# 
# f1 <- "./seurat_bam/scores_per_transcript.tab"
# f2 <- "./vipcca_bam/scores_per_transcript.tab"
# t1<-read.csv(f1, sep = "\t",quote = "\'")
# rownames(t1)<-markers
# t1<-t(t1[,4:16])
# ct1<-paste(rownames(t1),"S",sep = ".")
# t2<-read.csv(f2,sep = "\t",quote = "\'")
# rownames(t2)<-markers
# t2<-t(t2[,4:17])
# ct2<-paste(rownames(t2),"V",sep = ".")
# library(gtools)
# t<-smartbind(t1,t2,fill = 0)
# rownames(t)<-c(ct1,ct2)
# t<-t[order(rownames(t)),]
# t<-t[1:25,]
# m<-as.matrix(t(t))
# library("lattice")
# p3<-levelplot(m)


# lyz: ENSG00000090382:chr12:69742121-69748015
# gnly: ENSG00000115523: chr2:85912298-85925978
# ms4a1: ENSG00000156738:chr11:60223225-60238234
# bcl11b: ENSG00000127152: chr14:99635624-99737862
# cd8a: ENSG00000153563: chr2:87011729-87035520
# cd3e: ENSG00000198851: chr11:118175260-118186891
# cd4: ENSG00000010610: chr12:6896024-6929975
# cd14: ENSG00000170458: chr5:140011313-140013260


