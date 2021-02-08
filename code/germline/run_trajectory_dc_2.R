library(cowplot)
library(ggplot2)
library(Seurat)
library(monocle3)

rm(list=ls())
source("./code/run/evaluation_function.R")
# run trajectory on cvae results
result_path="./results/cell_align2/VIPCCA/"
glist <- list()
df_coef<-c()
cvae_list<-c()
cvae="cvae_10_128_12"
for(cvae in dir(result_path, pattern = "cvae")){
	curr_path=paste0(result_path,cvae,"/")
	print(curr_path)
	#h5ad_to_seurat(input_file = "output.h5ad", result_path = curr_path,
	#							 annotation_file="./data/cell_align/annotation_dc.txt")
	#datasets.integrated <- readRDS(paste0(curr_path,"integrated_data.rds"))
	#g<-plot_trajectory2(datasets.integrated, curr_path)
	#browser()

	cds <- readRDS(paste0(result_path,cvae,"/LPS & PAM_cds_pt_on_dr.rds"))
	traj <- unname(pseudotime(cds))
	colltime<-as.integer(substr(pData(cds)[["collection_time"]],1,1))
	coef<-cor.test(traj,colltime,method = "kendall")
	df_coef <- c(df_coef,coef$estimate)
	
	#g<-lapply(g, function(x) return(x+ggtitle(cvae)))
	#glist <- c(glist,g)
}
names(df_coef)<-dir(result_path, pattern = "cvae")

png(paste0(result_path,"all_umap.png"), width = 24*3, height =20*3, units = "in", res = 100)
plot_grid(plotlist = glist, ncol = 8, labels = "auto")
dev.off()

rm(list = ls())
source("./code/run/evaluation_function.R")

# run trajectory on uncorrected
result_path="./results/cell_align2/uncorrected/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory2(datasets.integrated, result_path, reduction = "pca")
result_path="./results/cell_align2/seurat/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory2(datasets.integrated, result_path, reduction = "pca")

# run trajectory on scanorama
result_path="./results/cell_align2/scanorama/"
h5ad_to_seurat(input_file = "output.h5ad", result_path = result_path,
							 annotation_file="./data/cell_align/annotation_dc.txt")
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory2(datasets.integrated, result_path, reduction = "scanorama")

# run trajectory on scalign
result_path="./results/cell_align2/scalign/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory2(datasets.integrated, result_path, reduction = "ALIGNED.CCA")

# run trajectory on MNN
result_path="./results/cell_align2/mnn/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory2(datasets.integrated, result_path, reduction = "pca")

# run trajectory on desc
result_path="./results/cell_align2/desc/"
h5ad_to_seurat(input_file = "output.h5ad", result_path = result_path,
							 annotation_file="./data/cell_align/annotation_dc.txt")
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
plot_trajectory2(datasets.integrated, result_path, reduction = "Embeded_z0.45")


##
rm(list = ls())
result_path="./results/cell_align2/uncorrected/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
dataset.pam<- subset(datasets.integrated, idents = "PAM")
#dataset.pam<-ScaleData(dataset.pam)
dataset.lps<- subset(datasets.integrated, idents = "LPS")
#dataset.lps<-ScaleData(dataset.lps)
cds.pam <- readRDS("./results/cell_align2/VIPCCA/s_cvae_10_128_12/PAM_cds_pt_on_dr.rds")
cds.lps <- readRDS("./results/cell_align2/VIPCCA/s_cvae_10_128_12/LPS_cds_pt_on_dr.rds")

dataset.pam@meta.data$order<-as.factor(rank(pseudotime(cds.pam)))
dataset.lps@meta.data$order<-as.factor(rank(pseudotime(cds.lps)))
Idents(dataset.pam)<-"order"
Idents(dataset.lps)<-"order"

library(patchwork)
e1<-DoHeatmap(dataset.pam, group.by = "collection_time",size = 50)+ggtitle("PAM")+NoLegend()+ylab("collection_time")
e2<-DoHeatmap(dataset.lps, group.by = "collection_time",size = 50)+ggtitle("LPS") +guides(color=FALSE,
                                                                                          fill = guide_colourbar(barwidth = 5, barheight = 15))+
  theme(legend.text = element_text(size = 40))
result_path="./results/cell_align2/"
result_dir<-c("uncorrected","desc","harmony","mnn","scalign","scanorama","seurat","VIPCCA/s_cvae_10_128_12")
method_names <- c("Before Integration","DESC","Harmony","MNN","scAlign","Scanorama", "Seurat V3", "VIPCCA")
ep<-e1+e2
ep<- ep&theme(title = element_text(size = 60))
df_coef<-data.frame()
for(i in 1:length(method_names)){
  res<-result_dir[i]
  print(paste0(result_path,res,"/PAM_cds_pt_on_dr.rds"))
  print(paste0(result_path,res,"/LPS_cds_pt_on_dr.rds"))
  cds.pam <- readRDS(paste0(result_path,res,"/PAM_cds_pt_on_dr.rds"))
  cds.lps <- readRDS(paste0(result_path,res,"/LPS_cds_pt_on_dr.rds"))
  dataset.pam@meta.data$order<-as.factor(rank(pseudotime(cds.pam)))
  dataset.lps@meta.data$order<-as.factor(rank(pseudotime(cds.lps))) 
  
  traj1<-unname(pseudotime(cds.pam))
  traj2<-unname(pseudotime(cds.lps))
  t1<-as.integer(substr(pData(cds.pam)[["collection_time"]],1,1))
  t2<-as.integer(substr(pData(cds.lps)[["collection_time"]],1,1))
  coef1<-cor.test(traj1,t1,method = "kendall")
  coef2<-cor.test(traj2,t2,method = "kendall")
  df_coef<-rbind(df_coef,c(coef1$estimate,coef2$estimate))
  
  e3<-DoHeatmap(dataset.pam,label = F,group.by = "order",draw.lines = FALSE)+ylab(method_names[i])+NoLegend()
  e4<-DoHeatmap(dataset.lps,label = F,group.by ="order", draw.lines = FALSE)+guides(color=FALSE)+NoLegend()
  ep<-ep+e3+e4
}
colnames(df_coef)<-c("pam.coef","lps.coef")
df_coef$method<-method_names
ep<-ep&theme(axis.title.y = element_text(size = 50))
png("./results/cell_align2/heatmap.png",width = 80, height = 60, units = "in", res = 100)
ep+plot_layout(ncol = 2)
dev.off()

df_coef_all<-data.frame(coef=c(df_coef$pam.coef,df_coef$lps.coef),
           meth=c(df_coef$method,df_coef$method),
           stm=rep(c("PAM","LPS"),each=8))
gcoef<-ggplot(df_coef_all,aes(x=meth,y=coef,group=stm))+
  geom_line(aes(color=stm),size=2)+ylab("Kendall coefficient")+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),
        text = element_text(size = 40))
png("./results/cell_align2/coef.png",width = 30, height = 15, units = "in",res = 100)
gcoef
dev.off()




# run trajectory on seurat
rm(list = ls())
result_path="./results/cell_align2/seurat/"
datasets.integrated <- readRDS(paste0(result_path,"integrated_data.rds"))
dataset.pam<- subset(datasets.integrated, idents = "PAM")
dataset.lps<- subset(datasets.integrated, idents = "LPS")
cds.lps<-readRDS("./results/cell_align2/seurat/LPS_cds_pt_on_dr.rds")
cds.pam<-readRDS("./results/cell_align2/seurat/PAM_cds_pt_on_dr.rds")
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


# run trajectory for comparison
rm(list = ls())
library(Seurat)
library(rlist)
source("./code/run/evaluation_function.R")
source("./code/cell_align2/cellAlign.R")
result_path="./results/cell_align2/"
tlist <- list()
glist <- list()
alist <- list()
result_dir<-c("uncorrected","desc","harmony","mnn","scalign","scanorama","seurat","VIPCCA/s_cvae_10_128_12")
reduction_names<-c("pca","Embeded_z0.45","harmony","pca","ALIGNED.CCA","scanorama", "pca","scxx")
method_names <- c("Before Integration","DESC","Harmony","MNN","scAlign","Scanorama", "Seurat V3", "VIPCCA")
i=0
df_corr_global<-data.frame()
df_corr_local<-data.frame()
for(res in result_dir){
	curr_path=paste0(result_path,res,"/")
	print(curr_path)
	i=i+1
	datasets.integrated <- readRDS(paste0(curr_path,"integrated_data.rds"))
	g1<-plot_trajectory(datasets.integrated, curr_path,
																t1="PAM",t2="LPS",method=method_names[i],reduction = reduction_names[i])

	g2<-plot_global_alignment3(curr_path,method = method_names[i])
	g3<-plot_local_alignment(curr_path,method = method_names[i])
	glist<-c(glist,g1,g2[1:2],g3[1])
	#browser()
	df_corr_global <- rbind(df_corr_global, g2[[3]])
	df_corr_local <- rbind(df_corr_local, g3[[2]])
}
# glist.copy<-glist
png("./results/cell_align2/corr_global.png",width = 15,height = 10,units = "in", res = 100)
gexp_global<-ggplot(df_corr_global, aes(x=corr,group=method,color=method))+geom_density(size=3)+theme_bw(base_size = 60)+
	theme(legend.position=c(0.4,0.6))+xlab("Correlation")
gexp_global
dev.off()

png("./results/cell_align2/corr_local.png",width = 15,height = 10,units = "in", res = 100)
gexp_local<-ggplot(df_corr_local, aes(x=corr,group=method,color=method))+geom_density(size=3)+theme_bw(base_size = 60)+
	#theme(legend.position=c(0.4,0.6))+
	xlab("Correlation")
gexp_local
dev.off()

glist<-lapply(glist, function(x) x+theme_bw(base_size = 40))
legend_traj<-cowplot::get_legend(glist[[1]]+guides(col=guide_legend("Collection time"))+
                                   theme(legend.text = element_text(size = 60),
                                                  legend.title = element_text(size = 60),
                                                  legend.direction = "horizontal"))
legend_stim<-cowplot::get_legend(glist[[4]]+theme(legend.title = element_blank(),
                                                  legend.text = element_text(size = 60),
                                                  legend.direction = "horizontal"))
glist<-lapply(glist, function(x) x+NoLegend())

library(rlist)
glist1<-glist[seq(1,40,5)]
glist2<-glist[seq(2,40,5)]
glist3<-glist[seq(4,40,5)]
glist4<-glist[seq(5,40,5)]

#glist3<-lapply(glist3, function(x)x+geom_point(size=10))

for(i in 1:8){
  glist3[[i]]<-glist3[[i]]+ggtitle(method_names[i])
  glist4[[i]]<-glist4[[i]]+ggtitle(method_names[i])
}

p1<-plot_grid(plotlist=glist1,ncol = 5)
p2<-plot_grid(plotlist=glist2,ncol = 5)
p1<-p1+ draw_label("PAM",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
p2<-p2+ draw_label("LPS",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size=60)+theme(plot.margin = margin(80, 0, 0, 0,unit = "pt"))
p3<-plot_grid(p1,p2,nrow=2,labels = c("a","b"),label_size = 60)
l1 <- plot_grid(p3,legend_traj,nrow = 2,rel_heights = c(8,1))

p4<-plot_grid(plotlist = glist3,ncol = 4)
p4<-p4+ draw_label("Global alignment of trajectories",fontface = 'bold',	x = 0.05,y=1.05,hjust = 0,size = 60)+theme(plot.margin = margin(100, 0, 0, 0,unit = "pt"))
p4 <- plot_grid(p4,legend_stim,nrow = 2,rel_heights = c(6,1))
p5<-gexp_global +theme(plot.margin = margin(100, 0, 100, 0,unit = "pt"))
l2 <- plot_grid(p4,p5,nrow = 1,rel_widths = c(4,2),labels = c("c","d"),label_size = 60)

p6<-plot_grid(l1,l2,nrow = 2,rel_heights = c(4,2))
png(paste0(result_path,"traj_inference.png"), width = 6*8, height =6*10, units = "in", res = 100)
p6
dev.off()

mode_pvalues<-compare_correlation(df_coef = df_corr_global)
mode_pvalues$logpvalue <- -log10(mode_pvalues$pvalue)
p7<-ggplot(mode_pvalues,aes(x=method,y=logpvalue))+
	geom_bar(stat = "identity",fill="orange")+
	ylab("-log(p-value)")+theme(axis.title.x = element_blank(),text = element_text(size=40))
png("./results/cell_align2/corr_gex_pvalue.png",width = 20,height = 20,units = "in",res = 100)
p7
dev.off()
