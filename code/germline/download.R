# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86146&format=file",
# 							destfile = "./data/germline/gse86146.rar")
# download.file("https://github.com/zorrodong/germcell/raw/master/germ_cell.RData",
# 								destfile = "./data/germline/germ_cell.RData")
# download.file("https://github.com/zorrodong/germcell/raw/master/TPM_and_UMI_counts.rar",
# 							destfile = "./data/germline/tpm_umi_counts.rar")
# download.file("https://zenodo.org/record/1443566/files/real/gold/germline-human-female-weeks_li.rds?download=1",
# 							destfile = "./data/germline/germline-female-li.rds")
# download.file("https://zenodo.org/record/1443566/files/real/gold/germline-human-male-weeks_li.rds?download=1",
# 							destfile = "./data/germline/germline-male-li.rds")
# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63818&format=file",
# 								destfile = "./data/germline/GSE63818_RAW.tar")
# download.file("https://zenodo.org/record/1443566/files/real/gold/germline-human-female_guo.rds?download=1",
# 							destfile = "./data/germline/germline-female-guo.rds")
# download.file("https://zenodo.org/record/1443566/files/real/gold/germline-human-male_guo.rds?download=1",
# 							destfile = "./data/germline/germline-male-guo.rds")

rm(list=ls())
fpkmToTpm <- function(fpkm)
{
	exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


i=1
for(fname in dir("./data/germline/GSE86146/",pattern="GSM",full.names = T)){
	x <- read.table(fname,header = T,row.names = 1)
	if(i>0){
		X<-x
		i=0
	}else{
		X<-cbind(X,x)
	}
}
X_li<-X
i=1
for(fname in dir("./data/germline/GSE63818_RAW/",pattern="GSM",full.names = T)){
	x <- read.table(fname,header = T,row.names = 1)
	if(i>0){
		X<-x
		i=0
	}else{
		X<-cbind(X,x)
	}
}
X_guo <- as.data.frame(apply(X,2,fpkmToTpm))

#colSums(X_guo)
#colSums(X_li)
# extract germ cells
cell.female.li<-readRDS("./data/germline/germline-female-li.rds")
cell.male.li <- readRDS("./data/germline/germline-male-li.rds")
cell.female.guo<-readRDS("./data/germline/germline-human-female_guo.rds")
cell.male.guo <- readRDS("./data/germline/germline-human-male_guo.rds")

cell.female.li.name<-cell.female.li$cell_ids
cell.male.li.name <- cell.male.li$cell_ids
cell.male.guo.name <- cell.male.guo$cell_info$title[grep(x=cell.male.guo$cell_info$title,pattern = "Soma",invert=T)]
cell.female.guo.name<-cell.female.guo$cell_info$title[grep(x=cell.female.guo$cell_info$title,pattern ="Soma",invert=T)]

feature.common<-intersect(rownames(X_guo),rownames(X_li))
feature.common.female <- intersect(cell.female.li$feature_info$feature_id,cell.female.guo$feature_info$feature_id)
feature.common.female <- intersect(feature.common,feature.common.female)
feature.common.male <- intersect(cell.male.li$feature_info$feature_id,cell.male.guo$feature_info$feature_id)
feature.common.male <- intersect(feature.common,feature.common.male)

X_female<-cbind(X_li[feature.common.female,cell.female.li.name],X_guo[feature.common.female,cell.female.guo.name])
X_male<-cbind(X_li[feature.common.male,cell.male.li.name],X_guo[feature.common.male,cell.male.guo.name])


week.female<-c(unlist(lapply(strsplit(cell.female.li.name,split = "_"),function(x) x[[2]])),
				unlist(lapply(strsplit(cell.female.guo.name,split = "_"),function(x) x[[3]])))
week.male<-c(unlist(lapply(strsplit(cell.male.li.name,split = "_"),function(x) x[[2]])),
							 unlist(lapply(strsplit(cell.male.guo.name,split = "_"),function(x) x[[3]])))
cell.female.info<-data.frame(batch=rep(c("Li","Guo"),times=c(length(cell.female.li.name),length(cell.female.guo.name))),rnames = colnames(X_female),week=week.female,row.names = "rnames")
cell.male.info <-data.frame(batch=rep(c("Li","Guo"),times=c(length(cell.male.li.name),length(cell.male.guo.name))),rnames=colnames(X_male),week=week.male,row.names = "rnames")


library(Seurat)
adata_female<-CreateSeuratObject(X_female,meta.data = cell.female.info)
adata_female@assays$RNA@data=log1p(adata_female@assays$RNA@counts)
adata_male<-CreateSeuratObject(X_male,meta.data = cell.male.info)
adata_male@assays$RNA@data=log1p(adata_male@assays$RNA@counts)

dir.create("./data/germline/female/",showWarnings = F)
dir.create("./data/germline/male/",showWarnings = F)
write.csv(X_female,"./data/germline/female/matrix.txt",quote = FALSE)
write.csv(X_male,"./data/germline/male/matrix.txt",quote = FALSE)
write.csv(cell.female.info,"./data/germline/female/annotation.txt",quote = FALSE)
write.csv(cell.male.info,"./data/germline/male/annotation.txt",quote = FALSE)
saveRDS(adata_female,"./data/germline/female/female.rds")
saveRDS(adata_male,"./data/germline/male/male.rds")

