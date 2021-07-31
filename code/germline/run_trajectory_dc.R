library(cowplot)
library(ggplot2)
library(Seurat)
library(monocle3)

rm(list=ls())
source("./code/run/evaluation_function.R")
# run trajectory on cvae results
result_path="./results/cell_align/VIPCCA/"
glist <- list()
for(cvae in dir(result_path, pattern = "cvae")){
	curr_path=paste0(result_path,cvae,"/")
	print(curr_path)
	h5ad_to_seurat(input_file = "output.h5ad", result_path = curr_path,
								 annotation_file="./data/cell_align/annotation_dc.txt")
	datasets.integrated <- readRDS(paste0(curr_path,"integrated_data.rds"))
	datasets.integrated <- run_umap_on_seurat(datasets.integrated, curr_path)
	g<-plot_trajectory(datasets.integrated, curr_path)
	#glist[[cvae]] <- g
	glist <- c(glist,g)
}
png(paste0(result_path,"all_umap.png"), width = 24*3, height =20*3, units = "in", res = 100)
plot_grid(plotlist = glist, ncol = 8, labels = "auto")
dev.off()

source("./code/run/evaluation_function.R")
rm(list = ls())
# run trajectory on uncorrected
result_path="./results/cell_align/uncorrected/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
dataset.pam<- subset(datasets.integrated, idents = "PAM")
#dataset.pam<-ScaleData(dataset.pam)
dataset.lps<- subset(datasets.integrated, idents = "LPS")
#dataset.lps<-ScaleData(dataset.lps)

#exp_pam<-
library(patchwork)
e1<-DoHeatmap(dataset.pam, group.by = "collection_time")+ggtitle("PAM")+NoLegend()
e2<-DoHeatmap(dataset.lps, group.by = "collection_time")+ggtitle("LPS")+guides(color=FALSE)
png("./results/cell_align/heatmap_expression_lag_2.png",width = 1500,height = 1000)
e1+e2
dev.off()
plot_trajectory(datasets.integrated, result_path, reduction = "pca")

# run trajectory on seurat
rm(list = ls())
result_path="./results/cell_align/seurat/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
dataset.pam<- subset(datasets.integrated, idents = "PAM")
dataset.lps<- subset(datasets.integrated, idents = "LPS")
cds.lps<-readRDS("./results/cell_align/seurat/LPS_cds_pt_on_dr.rds")
cds.pam<-readRDS("./results/cell_align/seurat/PAM_cds_pt_on_dr.rds")
df<-data.frame(
	medexpress=c(apply(dataset.lps@assays$integrated@scale.data,2,median),
							 apply(dataset.pam@assays$integrated@scale.data,2,median)),
	predtime=c(rescale(pseudotime(cds.lps)),rescale(pseudotime(cds.pam))),
	stim=rep(c("LPS","PAM"),times=c(dim(cds.lps)[2],dim(cds.pam)[2]))
	)
g<-ggplot(data = df, aes(x=predtime, y=medexpress)) + 
	geom_point(aes(colour=stim)) +
	ggtitle("tt")+ylab("Scaled expression")+xlab("Pseudotime")+
	theme(legend.title = element_blank(),
				legend.text = element_text(size=20))

plot_trajectory(datasets.integrated, result_path, reduction = "pca")

# run trajectory on scanorama
result_path="./results/cell_align/scanorama/"
h5ad_to_seurat(input_file = "output.h5ad", result_path = result_path,
							 annotation_file="./data/cell_align/annotation_dc.txt")
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory(datasets.integrated, result_path, reduction = "scanorama")

# run trajectory on scalign
result_path="./results/cell_align/scalign/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory(datasets.integrated, result_path, reduction = "ALIGNED.CCA")

# run trajectory on desc
result_path="./results/cell_align/desc/"
h5ad_to_seurat(input_file = "output.h5ad", result_path = result_path,
							 annotation_file="./data/cell_align/annotation_dc.txt")
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory(datasets.integrated, result_path, reduction = "Embeded_z0.45")

# run trajectory on liger
# result_path="./results/cell_align/liger/"
# h5ad_to_seurat(input_file = "output.h5ad", result_path = result_path,
# 							 annotation_file="./data/cell_align/annotation_dc.txt")
# datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
# plot_trajectory(datasets.integrated, result_path, reduction = "Embeded_z0.45")

# run trajectory for comparison
rm(list = ls())
library(Seurat)
library(rlist)
source("./code/run/evaluation_function.R")
source("./code/cell_align/cellAlign.R")
result_path="./results/cell_align/"
tlist <- list()
glist <- list()
alist <- list()
result_dir<-c("uncorrected","desc","harmony","scalign","scanorama","seurat","VIPCCA/cvae_10_64_16")
reduction_names<-c("pca","Embeded_z0.45","harmony","ALIGNED.CCA","scanorama", "pca","scxx")
method_names <- c("Before Integration","DESC","Harmony","scAlign","Scanorama", "Seurat V3", "VIPCCA")
i=0
# result_dir <- dir("./results/cell_align/VIPCCA/",pattern = "cvae")
datasets.integrated <- readRDS("./results/cell_align/uncorrected/integrated_data.rds")
dataset.pam<- subset(datasets.integrated, idents = "PAM")
dataset.lps<- subset(datasets.integrated, idents = "LPS")

medexpress<-c(apply(dataset.pam@assays$integrated@scale.data,2,median),
						 apply(dataset.lps@assays$integrated@scale.data,2,median))
mdexp1<-apply(dataset.pam@assays$integrated@scale.data,2,median)
mdexp2<-apply(dataset.lps@assays$integrated@scale.data,2,median)
for(res in result_dir){
	curr_path=paste0(result_path,res,"/")
	print(curr_path)
	i=i+1
	datasets.integrated <- readRDS(paste0(curr_path,"integrated_data.rds"))
	g1<-plot_trajectory(datasets.integrated, curr_path,
																t1="PAM",t2="LPS",method=method_names[i],reduction = reduction_names[i])
	legend_traj<-cowplot::get_legend(g1[[1]])
	#browser()
	g1[[1]]<-g1[[1]]+NoLegend()
	g1[[2]]<-g1[[2]]+NoLegend()
	
	
	# #
	#plot_alignment(curr_path)
	#
	cds1<-readRDS(paste0(curr_path,"PAM_cds_pt_on_dr.rds"))
	traj1<-unname(pseudotime(cds1))
	cds2<-readRDS(paste0(curr_path,"LPS_cds_pt_on_dr.rds"))
	traj2<-unname(pseudotime(cds2))
	t1<-as.integer(substr(pData(cds1)[["collection_time"]],1,1))
	t2<-as.integer(substr(pData(cds2)[["collection_time"]],1,1))
	
	set.seed(100)
	s1<-sample.int(length(traj1),500,replace = TRUE)
	s2<-sample.int(length(traj2),500,replace = TRUE)
	traj1.sample <- traj1[s1]
	mdexp1.sample <- mdexp1[s1]
	traj2.sample <- traj2[s2]
	mdexp2.sample <- mdexp2[s2]
	
	cor.res<-cor.test(mdexp1.sample[order(traj1.sample)],
										mdexp2.sample[order(traj2.sample)])
	
	tlist[[res]] <- c(cor.res$estimate,cor.res$p.value)
	
	df<-data.frame(
		medexpress=medexpress,
		predtime=c(rescale(traj1),rescale(traj2)),
		stim=rep(c("PAM","LPS"),times=c(length(traj1),length(traj2)))
	)
	g2<-ggplot(data = df, aes(x=predtime, y=medexpress, colour=stim)) + 
		geom_point(aes(colour=stim)) + geom_smooth(se = T)+
		ggtitle(method_names[i])+ylab("Scaled expression")+xlab("Pseudotime")+
		theme(legend.title = element_blank(),
					legend.text = element_text(size=20),
					legend.position = c(0.4,0.8))
	legend_stim<-cowplot::get_legend(g2)
	g2<-g2+NoLegend()
	glist<-c(glist,g1)
	glist<-list.append(glist,g2)
	#model<-
	#model <- lm(sales ~ youtube, data = marketing)
	
}

df <- as.data.frame(t(as.data.frame(tlist)))
colnames(df) <- c("Pearson","pvalue")
df$pvalue<- -log10(df$pvalue)
g<-ggplot(data = df, aes(x=pvalue, y=Pearson)) + 
	geom_point(aes(colour=method_names)) +
	ggtitle("Kendall correlation coefficient")+
	theme(legend.title = element_blank(),
				legend.text = element_text(size=20))
legend1 <- cowplot::get_legend(g)
#g<-g+NoLegend()

#galign_raw <- image_read("./results/cell_align/VIPCCA/cvae_10_64_16/global_alignment.png")
library(rlist)
glist1<-glist[seq(1,18,3)]
glist2<-glist[seq(2,18,3)]
glist3<-glist[seq(3,19,3)]

#glist2<-list.append(glist,g,legend1)
#alist<-alist[c(seq(1,12,2),seq(2,12,2))]

title1 <- ggdraw() + draw_label("PAM",fontface = 'bold',	x = 0,hjust = 0) +theme(plot.margin = margin(0, 0, 0, 7))
title2 <- ggdraw() + draw_label("LPS",fontface = 'bold',	x = 0,hjust = 0) +theme(plot.margin = margin(0, 0, 0, 7))


p1<-plot_grid(plotlist=glist1,ncol = 2,labels="auto",label_size = 30)
p2<-plot_grid(plotlist=glist2,ncol = 2,labels=c("g","h","i","j","k","l"),label_size = 30)

p3<-plot_grid(title1,title2,p1,p2,ncol=2,
							rel_heights = c(0.1,1))
p4<-plot_grid(plotlist = glist3,ncol = 4,labels = c("m","n","o","p","q","r"),label_size = 30)
p5<-plot_grid(p3,legend_traj,p4,legend_stim,ncol = 2,rel_heights = c(3,2),rel_widths = c(6,0.4))
png(paste0(result_path,"traj_inference.png"), width = 4*5, height =4*5, units = "in", res = 300)
p5
dev.off()
