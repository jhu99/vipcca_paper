rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)

#### Download 
urls<-c(
	"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5FcDNA%2Ebarcodes%2Etsv%2Egz",
	"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5FcDNA%2Ecounts%2Emtx%2Egz",
	"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5FcDNA%2Egenes%2Etsv%2Egz",
	"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5Fchromatin%2Ebarcodes%2Etsv%2Egz",
	"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5Fchromatin%2Ecounts%2Emtx%2Egz",
	"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126074&format=file&file=GSE126074%5FP0%5FBrainCortex%5FSNAREseq%5Fchromatin%2Epeaks%2Etsv%2Egz"
	)
fnames<-c(
	"GSE126074_P0_BrainCortex_SNAREseq_cDNA/barcodes.tsv.gz",
	"GSE126074_P0_BrainCortex_SNAREseq_cDNA/matrix.mtx.gz",
	"GSE126074_P0_BrainCortex_SNAREseq_cDNA/features.tsv.gz",
	"GSE126074_P0_BrainCortex_SNAREseq_chromatin/barcodes.tsv.gz",
	"GSE126074_P0_BrainCortex_SNAREseq_chromatin/matrix.mtx.gz",
	"GSE126074_P0_BrainCortex_SNAREseq_chromatin/features.tsv.gz"
)
# for(i in 1:6){
# 	download.file(url = urls[i],destfile = paste0("./data/atacseq/snareseq/",fnames[i]))
# }
# mm10gtf="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz"
# download.file(mm10gtf,"./data/atacseq/snareseq/mm10_refgene.gtf.gz")

#### create gene activity matrix
peak<-Read10X("./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_chromatin/",gene.column = 1)
# geneActMat<-CreateGeneActivityMatrix(peak.matrix = peak,
# 																		 annotation.file = "./data/atacseq/snareseq/Mus_musculus.GRCm38.100.gtf",
# 																		 seq.levels = c(1:19,"X","Y")
# 																		 )
# saveRDS(geneActMat,"./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_chromatin/geneActMat.rds")
geneActMat<-readRDS("./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_chromatin/geneActMat.rds")
ExpressionMat<-Read10X("./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_cDNA/",gene.column = 1)
bcmeta<-readRDS("./data/atacseq/snareseq/P0BrainCortex_SNAREseq_metadata_full.rds")
geneset1<-rownames(geneActMat)
geneset2<-rownames(ExpressionMat)
genelist<-read.csv("./data/atacseq/snareseq/12A/genelist_python.txt",header = FALSE,col.names = "hvg")$hvg
commongenes<-intersect(geneset1,geneset2)
rownames(bcmeta) <- bcmeta$Barcode
## Seurat object setup
bcortex.atac<-CreateSeuratObject(counts = peak, meta.data = bcmeta, assay = "ATAC", project = "bcortex_atac")
bcortex.atac[["ACTIVITY"]]<-CreateAssayObject(counts = geneActMat[commongenes,])
bcortex.atac$tech <- "atac"
Idents(bcortex.atac)<-"Batch"
bcortex.atac<-subset(bcortex.atac, idents = "12A" )
library(Matrix)
writeMM(bcortex.atac@assays$ACTIVITY@counts,"./data/atacseq/snareseq/12A/GeneActivity_chromatin/matrix.mtx")
write.table(colnames(bcortex.atac),"./data/atacseq/snareseq/12A/GeneActivity_chromatin/barcodes.tsv",
						row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(cbind(commongenes,commongenes),"./data/atacseq/snareseq/12A/GeneActivity_chromatin/genes.tsv",sep="\t",
						row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bcortex.atac@meta.data,"./data/atacseq/snareseq/annotation.csv")
bcortex.rna<-CreateSeuratObject(counts = ExpressionMat[commongenes,], meta.data = bcmeta, assay = "RNA", project = "bcortex_rna")
bcortex.rna$tech <-"scRNAseq"
Idents(bcortex.rna)<-"Batch"
bcortex.rna<-subset(bcortex.rna, idents = "12A" )
writeMM(bcortex.rna@assays$RNA@counts,"./data/atacseq/snareseq/12A/ExpressionMatrix_cDNA/matrix.mtx")
write.table(colnames(bcortex.rna),"./data/atacseq/snareseq/12A/ExpressionMatrix_cDNA/barcodes.tsv",
						row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(cbind(commongenes,commongenes),"./data/atacseq/snareseq/12A/ExpressionMatrix_cDNA/genes.tsv",sep = "\t",
						row.names = FALSE, col.names = FALSE, quote = FALSE)

## preprocessing ATAC and geneActivity
DefaultAssay(bcortex.atac) <- "ACTIVITY"
bcortex.atac <- NormalizeData(bcortex.atac,scale.factor = 1e6)
VariableFeatures(bcortex.atac)<-genelist
bcortex.atac <- FindVariableFeatures(bcortex.atac,selection.method = "vst", nfeatures = 2000)
bcortex.atac <- ScaleData(bcortex.atac,features = rownames(bcortex.atac))
DefaultAssay(bcortex.atac) <- "ATAC"
VariableFeatures(bcortex.atac) <- names(which(Matrix::rowSums(bcortex.atac) > 100))
bcortex.atac <- RunLSI(bcortex.atac, n = 50, scale.max = NULL)
bcortex.atac <- RunUMAP(bcortex.atac, reduction = "lsi", dims = 1:50)
p1 <- DimPlot(bcortex.atac, reduction = "umap",group.by = "IdentChar") + NoLegend() + ggtitle("ATAC")
saveRDS(bcortex.atac,"./data/atacseq/snareseq/12A/GeneActivity_chromatin/atac_sobj.rds")
# Integrating rna from replicates
# obj.list<-SplitObject(bcortex.rna,split.by = "Batch")
# for (i in 1:length(obj.list)) {
# 	obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE, scale.factor = 1e4)
# 	obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst",
# 	 																					 nfeatures = 2000, verbose = FALSE)
# }
# obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:16)
# obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:16)
# DefaultAssay(obj.integrated) <- "integrated"
# obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
# obj.integrated <- RunPCA(obj.integrated, npcs = 16, verbose = FALSE)
# obj.integrated <- RunUMAP(obj.integrated,reduction = "pca",dims = 1:16)
# p1<-DimPlot(obj.integrated,reduction = "umap", group.by = "Batch")
# p2<-DimPlot(obj.integrated,reduction = "umap", group.by = "IdentChar")
# p1+p2

bcortex.rna <- NormalizeData(bcortex.rna,scale.factor = 1e6)
VariableFeatures(bcortex.rna)<-genelist
bcortex.rna <- FindVariableFeatures(bcortex.rna,selection.method = "vst", nfeatures = 2000)
bcortex.rna <- ScaleData(bcortex.rna, features=rownames(bcortex.rna))
bcortex.rna <- RunPCA(bcortex.rna, npcs = 30)

bcortex.rna <- RunUMAP(bcortex.rna,reduction = "pca", dims = 1:30)
features <- c("Slc1a3",#"Mik67","Ararb2",
							"Gadd45g","Unc5d","Cntn2","Rbfox1","Tenm3","Foxp1","Crmp1",
							"Galntl6", "Epha6", "Sox5", "Nxph1", "Reln","Pdgfra","Cldn5","Cald1","Apbb1ip"
							)
#p3 <- DimPlot(bcortex.rna, reduction = "umap", group.by = "Batch") 
p4 <- DimPlot(bcortex.rna, reduction = "umap", group.by = "IdentChar")+ggtitle("RNA")
p1+p4

p5 <- FeaturePlot(bcortex.rna, features = features, ncol = 6)
saveRDS(bcortex.rna,"./data/atacseq/snareseq/12A/ExpressionMatrix_cDNA/rna_sobj.rds")
png("./results/atacseq/12A/paired_plot.png",width = 16, height = 8,units = "in", res = 150)
p1 +p4
dev.off()
png("./results/atacseq/12A/feature_plot.png",width = 32, height = 16,units = "in", res = 150)
p5
dev.off()
table(bcortex.rna$Batch,bcortex.rna$IdentChar)
transfer.anchors <- FindTransferAnchors(reference = bcortex.rna, query = bcortex.atac, features = VariableFeatures(object = bcortex.rna), 
																				reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = bcortex.rna$IdentChar, 
																		 weight.reduction = bcortex.atac[["lsi"]])
bcortex.atac <- AddMetaData(bcortex.atac, metadata = celltype.predictions)
png("./results/atacseq/12A/seurat_prediction_hist.png")
hist(bcortex.atac$prediction.score.max)
dev.off()
celltype.predictions$celltype=bcortex.atac$IdentChar
#celltype.predictions<-celltype.predictions[celltype.predictions$prediction.score.max>0.5,]
indicator=(celltype.predictions$celltype==celltype.predictions$predicted.id)
sum(indicator)/length(indicator)
saveRDS(celltype.predictions,"./results/atacseq/12A/seurat.celltype.predictions.rds")
# bcortex.rna$seurat_clusters
# f<-read.table("./data/atacseq/snareseq/GSE126074_P0_BrainCortex_SNAREseq_cDNA/features.tsv.gz")
rm(list = ls())
library(Seurat)
library(patchwork)
ref="RNA"
query="ATAC"
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
fd="./results/atacseq/12A/scxx_cvae/scmb_cvae_2_64_16"
for(fd in dir("./results/atacseq/12A/scxx_cvae",full.names = T)){
	adata<-ReadH5AD(paste(fd,"output.h5ad",sep = "/"))
	ann<-read.csv("./data/atacseq/snareseq/12A/annotation_vipcca.csv",row.names = 1)
	adata<-AddMetaData(adata,metadata = ann)
	prediction<-celltypePrediction(adata,ref="RNA",query="ATAC")
	prediction_filtered<-prediction[prediction$score>0.0,]
	indicators<-(as.character(prediction_filtered$prediction)==as.character(prediction_filtered$celltype))
	print(fd)
	print(sum(indicators)/length(indicators))
	# calculate a prediction score
	# visulization
	# adata<-RunUMAP(adata, reduction = "scxx", dims = 1:16, min.dist = 0.1,)
	# p6<-DimPlot(adata, reduction = "umap", group.by = "X_batch")
	# p7<-DimPlot(adata, reduction = "umap", group.by = "IdentChar")
	# png(paste(fd,"paired_plot.png",sep="/"),width = 16, height = 8, units = "in", res = 150)
	# print(p6+p7)
	# dev.off()
}

#x<-read.csv("./results/1M_mouse_brain/cpu_mem_profile_scanorama.txt")
