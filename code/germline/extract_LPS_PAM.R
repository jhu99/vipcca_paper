# Extract LPS and PAM
rm(list = ls())
x1<-read.csv("./data/cell_align/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt",sep = "\t")
xpam<-read.csv("./data/cell_align/pam.txt",header = T,row.names = 1)
xlps<-read.csv("./data/cell_align/lps.txt",header = T,row.names = 1)
x1.selected.pam <- x1[,colnames(xpam)]
x1.selected.lps <- x1[,colnames(xlps)]
write.csv(x1.selected.pam,"./data/cell_align/GSE48968_allgenes/pam.txt")
write.csv(x1.selected.lps,"./data/cell_align/GSE48968_allgenes/lps.txt")
# download gene modules
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature13437/MediaObjects/41586_2014_BFnature13437_MOESM107_ESM.xls",
							destfile = "./data/cell_align/GSE48968_allgenes/gene_modules.xls")
genemodules<-readxl::read_excel("./data/cell_align/GSE48968_allgenes/gene_modules.xls",
																col_types = c("text","text","text"))
genelist<-genemodules$GENE
genelist[185]="MARCH5"
genemodules[185,1]="MARCH5"
intersect(genelist,rownames(x1))
##run scanpy.pp.highly_variable_genes
## code in prepare_test_data.py
## xxx


## 

genelist
x2<- x1[genelist,colnames(xpam)]
x3<- x1[genelist,colnames(xlps)]
sum(is.na(x2))
write.csv(x2,file = "./data/cell_align/GSE48968_allgenes/pam_gm813.txt")
write.csv(x3,file = "./data/cell_align/GSE48968_allgenes/lps_gm813.txt")

# output gene module.
gmid <- genemodules[genemodules$CLUSTER=="Id",]
gmiiic <- genemodules[genemodules$CLUSTER=="IIIc",]
x2<- x1[gmid$GENE,colnames(xpam)]
x3<- x1[gmid$GENE,colnames(xlps)]
write.csv(x2,file = "./data/cell_align/GSE48968_allgenes/pam_gmid.txt")
write.csv(x3,file = "./data/cell_align/GSE48968_allgenes/lps_gmid.txt")
x2<- x1[gmiiic$GENE,colnames(xpam)]
x3<- x1[gmiiic$GENE,colnames(xlps)]
write.csv(x2,file = "./data/cell_align/GSE48968_allgenes/pam_gmiiic.txt")
write.csv(x3,file = "./data/cell_align/GSE48968_allgenes/lps_gmiiic.txt")


# x=seq(0,10,0.01)
# a=8
# y=1/(1+x)^0.5+1/(1+a)^0.5+(a*x/(a*x+8))^0.5
# ggplot(data.frame(x=x,y=y),aes(x=x,y=y))+
# 	geom_line()

