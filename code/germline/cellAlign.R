library(cellAlign)
library(monocle3)
# load data and trajectory

plot_alignment<-function(result_path="./results/cell_align/VIPCCA/cvae_10_64_16/"){
	data(expGlobalLPS)
	data(expGlobalPAM)
	
	cds1<-readRDS(paste0(result_path,"LPS_cds_pt_on_dr.rds"))
	cds2<-readRDS(paste0(result_path,"PAM_cds_pt_on_dr.rds"))
	trajLPS <- pseudotime(cds1)
	trajPAM <- pseudotime(cds2)
	trajLPS<-trajLPS/max(trajLPS)
	trajPAM<-trajPAM/max(trajPAM)
	# interplolation
	
	numPts = 200
	interGlobalLPS = cellAlign::interWeights(expDataBatch = expGlobalLPS, trajCond = trajLPS,
																					 winSz = 0.1, numPts = numPts)
	interGlobalPAM = cellAlign::interWeights(expDataBatch = expGlobalPAM, trajCond = trajPAM,
																					 winSz = 0.1, numPts = numPts)
	
	require(ggplot2)
	require(reshape2)
	require(pheatmap)
	
	sharedMarkers = intersect(rownames(expGlobalLPS), rownames(expGlobalPAM))
	whichgene=sharedMarkers[1]
	selectedLPS<-interGlobalLPS$interpolatedVals[whichgene,]
	selectedPAM<-interGlobalPAM$interpolatedVals[whichgene,]
	
	dfLPSi = data.frame(traj = interGlobalLPS$traj, value=(selectedLPS), error=interGlobalLPS$error[whichgene,])
	dfLPS = data.frame(traj = trajLPS, t(expGlobalLPS[whichgene,]))
	dfPAMi = data.frame(traj = interGlobalPAM$traj, value=(selectedPAM), error=interGlobalPAM$error[whichgene,])
	dfPAM = data.frame(traj = trajPAM, t(expGlobalPAM[whichgene,]))
	dfLPSM = melt(dfLPS, id.vars = 'traj')
	dfPAMM = melt(dfPAM, id.vars = 'traj')
	#plot of an example gene and its interpolation with error bars
	#ggplot(dfLPSi, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=dfLPSM, aes(x=traj,y=value)) + ggtitle(whichgene) 
	
	#scale the interpolated data (Recommended):
	interScaledGlobalLPS = cellAlign::scaleInterpolate(interGlobalLPS)
	interScaledGlobalPAM = cellAlign::scaleInterpolate(interGlobalPAM)
	#calculate dissimilar matrix
	# A=calcDistMat(interScaledGlobalPAM$scaledData[,1:10],interScaledGlobalLPS$scaledData[,1:10], dist.method = 'Euclidean')
	# pheatmap(A, cluster_cols = F, cluster_rows=F, main = "LPS vs PAM distances, 1st 10 points",
	# 				 show_rownames = F, show_colnames = F,display_numbers = TRUE)
	#global alignment
	alignment = globalAlign(interScaledGlobalPAM$scaledData, interScaledGlobalLPS$scaledData,
													scores = list(query = interScaledGlobalPAM$traj, 
																				ref = interScaledGlobalLPS$traj),
													sigCalc = F, numPerm = 20)
  #browser()
	
	png(paste0(result_path,"global_alignment_test.png"), width = 800, height = 600)
	plotAlign(alignment)
	dev.off()
	
	# mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalPAM$traj, realTrajQuery = trajPAM,
	# 														intTrajRef = interScaledGlobalLPS$traj, realTrajRef = trajLPS)
	# png(paste0(result_path,"global_mapping.png"), width = 800, height = 600)
	# g2<-plotMapping(mapping)
	# print(g2)
	# dev.off()
}

plot_alignment2<-function(result_path="./results/cell_align/VIPCCA/cvae_10_64_16/"){
  data(expGlobalLPS)
  data(expGlobalPAM)
  cds1<-readRDS(paste0(result_path,"LPS_cds_pt_on_dr.rds"))
  cds2<-readRDS(paste0(result_path,"PAM_cds_pt_on_dr.rds"))
  trajLPS <- pseudotime(cds1)
  trajPAM <- pseudotime(cds2)
  trajLPS<-trajLPS/max(trajLPS)
  trajPAM<-trajPAM/max(trajPAM)
  # interplolation
  #browser()
  numPts = 200
  interGlobalLPS = cellAlign::interWeights(expDataBatch = expGlobalLPS, trajCond = trajLPS,
                                           winSz = 0.1, numPts = numPts)
  interGlobalPAM = cellAlign::interWeights(expDataBatch = expGlobalPAM, trajCond = trajPAM,
                                           winSz = 0.1, numPts = numPts)
  #scale the interpolated data (Recommended):
  interScaledGlobalLPS = cellAlign::scaleInterpolate(interGlobalLPS)
  interScaledGlobalPAM = cellAlign::scaleInterpolate(interGlobalPAM)

  alignment = globalAlign(interScaledGlobalPAM$scaledData, interScaledGlobalLPS$scaledData,
                          scores = list(query = interScaledGlobalPAM$traj, 
                                        ref = interScaledGlobalLPS$traj),
                          sigCalc = F, numPerm = 20)
  
  mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalPAM$traj, realTrajQuery = trajPAM,
  														intTrajRef = interScaledGlobalLPS$traj, realTrajRef = trajLPS)
  png(paste0(result_path,"global_mapping.png"), width = 800, height = 600)
  g1<-plotMapping(mapping)
  print(g1)
  dev.off()

  mdsLPS<-data.frame(mdexp=colMedians(interScaledGlobalLPS$scaledData),
                     traj=(1:200)/200,cond="LPS")
  mdsPAM<-data.frame(mdexp=colMedians(interScaledGlobalPAM$scaledData),
                     traj=(1:200)/200,cond="PAM")
  
  mdsALL<-c()
  # browser()
  numConn<-length(alignment$align[[1]]$index1)
  for(i in 1:numConn){
    ppam <- alignment$align[[1]]$index1[i]
    plps <- alignment$align[[1]]$index2[i]
    com <- rbind(mdsLPS[plps,],mdsPAM[ppam,])
    com$paired <- i
    mdsALL<-rbind(mdsALL,com)
  }
  mdsALL<-mdsALL[c(seq(1,2*numConn,2),seq(2,2*numConn,2)),]
  s<-seq(1,numConn,4)
  mdsALL<-mdsALL[c(s,s+numConn),]
  g2<-ggplot(mdsALL,aes(x=traj,y=mdexp,group=cond,color=cond))+
    geom_point(size=6)+geom_line(aes(group=paired),size=2,color="grey")+theme_bw(base_size = 40)+
    xlab("Scaled pseudotime")+ylab("Scaled expression")+
    theme(legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.text = element_text(size=60),
          legend.position = c(0.4,0.8))
  
  return(list(g1,g2))
}

plot_global_alignment3<-function(result_path="./results/cell_align2/VIPCCA/cvae_10_128_20/",method="VIPCCA"){
	xlps <- read.csv("./data/cell_align/GSE48968_allgenes/lps_gmid.txt",header = T,row.names = 1)
	xpam <- read.csv("./data/cell_align/GSE48968_allgenes/pam_gmid.txt",header = T,row.names = 1)
	xlps <- log1p(xlps)
	xpam <- log1p(xpam)
	
	cds1<-readRDS(paste0(result_path,"LPS_cds_pt_on_dr.rds"))
	cds2<-readRDS(paste0(result_path,"PAM_cds_pt_on_dr.rds"))
	trajLPS <- pseudotime(cds1)
	trajPAM <- pseudotime(cds2)
	trajLPS<-trajLPS/max(trajLPS)
	trajPAM<-trajPAM/max(trajPAM)
	# interplolation
	#browser()
	numPts = 200
	interGlobalLPS = cellAlign::interWeights(expDataBatch = xlps, trajCond = trajLPS,
																					 winSz = 0.1, numPts = numPts)
	interGlobalPAM = cellAlign::interWeights(expDataBatch = xpam, trajCond = trajPAM,
																					 winSz = 0.1, numPts = numPts)
	#scale the interpolated data (Recommended):
	interScaledGlobalLPS = cellAlign::scaleInterpolate(interGlobalLPS)
	interScaledGlobalPAM = cellAlign::scaleInterpolate(interGlobalPAM)
	
	alignment = globalAlign(interScaledGlobalPAM$scaledData, interScaledGlobalLPS$scaledData,
													scores = list(query = interScaledGlobalPAM$traj, 
																				ref = interScaledGlobalLPS$traj),
													sigCalc = F, numPerm = 20)
	
	mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalPAM$traj, realTrajQuery = trajPAM,
															intTrajRef = interScaledGlobalLPS$traj, realTrajRef = trajLPS)
	png(paste0(result_path,"global_mapping.png"), width = 800, height = 600)
	g1<-plotMapping(mapping)
	print(g1)
	dev.off()
	
	mdsLPS<-data.frame(mdexp=colMedians(interScaledGlobalLPS$scaledData),
										 traj=(1:200)/200,cond="LPS")
	mdsPAM<-data.frame(mdexp=colMedians(interScaledGlobalPAM$scaledData),
										 traj=(1:200)/200,cond="PAM")
	
	mdsALL<-c()
	# browser()
	numConn<-length(alignment$align[[1]]$index1)
	for(i in 1:numConn){
		ppam <- alignment$align[[1]]$index1[i]
		plps <- alignment$align[[1]]$index2[i]
		com <- rbind(mdsLPS[plps,],mdsPAM[ppam,])
		com$paired <- i
		mdsALL<-rbind(mdsALL,com)
		
	}
	# It calcualates correlation of gene expression patterns from aligned trajectories here.
	xpam_aligned<-interScaledGlobalPAM$scaledData[,alignment$align[[1]]$index1]
	xlps_aligned<-interScaledGlobalLPS$scaledData[,alignment$align[[1]]$index2]
	library(rlist)
	coef_genes<-list()
	for(i in 1:dim(xpam_aligned)[1]){
		coef<-cor.test(xpam_aligned[i,],xlps_aligned[i,],method = "pearson")
		coef_genes<-list.append(coef_genes,coef)
	}
	coef_genes<-unlist(lapply(coef_genes,function(x) return(x$estimate)))
	coef_genes<-unname(unlist(coef_genes))
	df_coef<-data.frame(corr=coef_genes,method=rep(method,length(coef_genes)))
	g<-ggplot(df_coef,aes(x=corr))+geom_density()
	
	mdsALL<-mdsALL[c(seq(1,2*numConn,2),seq(2,2*numConn,2)),]
	s<-seq(1,numConn,4)
	mdsALL<-mdsALL[c(s,s+numConn),]
	g2<-ggplot(mdsALL,aes(x=traj,y=mdexp,group=cond,color=cond))+
		geom_point(size=6)+geom_line(aes(group=paired),size=2,color="grey")+theme_bw(base_size = 40)+
		xlab("Scaled pseudotime")+ylab("Scaled expression")+
		theme(legend.title = element_blank(),
					legend.direction = "horizontal",
					legend.text = element_text(size=60),
					legend.position = c(0.4,0.8))
	
	return(list(g1,g2,df_coef))
}

plot_global_alignment4<-function(result_path="./results/cell_align2/VIPCCA/cvae_10_128_12/",method="VIPCCA"){
	xlps <- read.csv("./data/cell_align/GSE48968_allgenes/lps_gmid.txt",header = T,row.names = 1)
	xpam <- read.csv("./data/cell_align/GSE48968_allgenes/pam_gmid.txt",header = T,row.names = 1)
	xlps <- log1p(xlps)
	xpam <- log1p(xpam)
	
	cds <- readRDS(paste0(result_path,"LPS & PAM_cds_pt_on_dr.rds"))
	cds1<-cds[,pData(cds)$condition=="LPS"]
	cds2<-cds[,pData(cds)$condition=="PAM"]
	trajLPS <- pseudotime(cds1)
	trajPAM <- pseudotime(cds2)
	trajLPS<-trajLPS/max(c(trajLPS,trajPAM))
	trajPAM<-trajPAM/max(c(trajLPS,trajPAM))
	# interplolation
	#browser()
	numPts = 200
	interGlobalLPS = cellAlign::interWeights(expDataBatch = xlps, trajCond = trajLPS,
																					 winSz = 0.1, numPts = numPts)
	interGlobalPAM = cellAlign::interWeights(expDataBatch = xpam, trajCond = trajPAM,
																					 winSz = 0.1, numPts = numPts)
	#scale the interpolated data (Recommended):
	interScaledGlobalLPS = cellAlign::scaleInterpolate(interGlobalLPS)
	interScaledGlobalPAM = cellAlign::scaleInterpolate(interGlobalPAM)
	
	alignment = globalAlign(interScaledGlobalPAM$scaledData, interScaledGlobalLPS$scaledData,
													scores = list(query = interScaledGlobalPAM$traj, 
																				ref = interScaledGlobalLPS$traj),
													sigCalc = F, numPerm = 20)
	
	mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalPAM$traj, realTrajQuery = trajPAM,
															intTrajRef = interScaledGlobalLPS$traj, realTrajRef = trajLPS)
	png(paste0(result_path,"global_mapping.png"), width = 800, height = 600)
	g1<-plotMapping(mapping)
	print(g1)
	dev.off()
	
	mdsLPS<-data.frame(mdexp=colMedians(interScaledGlobalLPS$scaledData),
										 traj=(1:200)/200,cond="LPS")
	mdsPAM<-data.frame(mdexp=colMedians(interScaledGlobalPAM$scaledData),
										 traj=(1:200)/200,cond="PAM")
	
	mdsALL<-c()
	# browser()
	numConn<-length(alignment$align[[1]]$index1)
	for(i in 1:numConn){
		ppam <- alignment$align[[1]]$index1[i]
		plps <- alignment$align[[1]]$index2[i]
		com <- rbind(mdsLPS[plps,],mdsPAM[ppam,])
		com$paired <- i
		mdsALL<-rbind(mdsALL,com)
		
	}
	# It calcualates correlation of gene expression patterns from aligned trajectories here.
	xpam_aligned<-interScaledGlobalPAM$scaledData[,alignment$align[[1]]$index1]
	xlps_aligned<-interScaledGlobalLPS$scaledData[,alignment$align[[1]]$index2]
	library(rlist)
	coef_genes<-list()
	for(i in 1:dim(xpam_aligned)[1]){
		coef<-cor.test(xpam_aligned[i,],xlps_aligned[i,],method = "pearson")
		coef_genes<-list.append(coef_genes,coef)
	}
	coef_genes<-unlist(lapply(coef_genes,function(x) return(x$estimate)))
	coef_genes<-unname(unlist(coef_genes))
	df_coef<-data.frame(corr=coef_genes,method=rep(method,length(coef_genes)))
	g<-ggplot(df_coef,aes(x=corr))+geom_density()
	
	mdsALL<-mdsALL[c(seq(1,2*numConn,2),seq(2,2*numConn,2)),]
	s<-seq(1,numConn,4)
	mdsALL<-mdsALL[c(s,s+numConn),]
	g2<-ggplot(mdsALL,aes(x=traj,y=mdexp,group=cond,color=cond))+
		geom_point(size=6)+geom_line(aes(group=paired),size=2,color="grey")+theme_bw(base_size = 40)+
		xlab("Scaled pseudotime")+ylab("Scaled expression")+
		theme(legend.title = element_blank(),
					legend.direction = "horizontal",
					legend.text = element_text(size=60),
					legend.position = c(0.4,0.8))
	
	return(list(g1,g2,df_coef))
}

plot_local_alignment<-function(result_path="./results/cell_align2/VIPCCA/s_cvae_10_128_12/",method="VIPCCA"){
	xlps <- read.csv("./data/cell_align/GSE48968_allgenes/lps_gmiiic.txt",header = T,row.names = 1)
	xpam <- read.csv("./data/cell_align/GSE48968_allgenes/pam_gmiiic.txt",header = T,row.names = 1)
	xlps <- log1p(xlps)
	xpam <- log1p(xpam)
	
	cds <- readRDS(paste0(result_path,"LPS & PAM_cds_pt_on_dr.rds"))
	cds1<-cds[,pData(cds)$condition=="LPS"]
	cds2<-cds[,pData(cds)$condition=="PAM"]
	trajLPS <- pseudotime(cds1)
	trajPAM <- pseudotime(cds2)
	trajLPS<-trajLPS/max(c(trajLPS,trajPAM))
	trajPAM<-trajPAM/max(c(trajLPS,trajPAM))
	# interplolation
	#browser()
	
	numPts=200
	Thresh=0.3
	interLocalLPS = interWeights(expDataBatch = xlps, trajCond = trajLPS, winSz = 0.1, numPts = numPts)
	interLocalPAM = interWeights(expDataBatch = xpam, trajCond = trajPAM, winSz = 0.1, numPts = numPts)
	interScaledLocalLPS = cellAlign::scaleInterpolate(interLocalLPS)
	interScaledLocalPAM = cellAlign::scaleInterpolate(interLocalPAM)
	
	A=calcDistMat(interScaledLocalPAM$scaledData,interScaledLocalLPS$scaledData, dist.method = 'Euclidean')
	A[A > 10*Thresh] <- max(A)
	alignment = localAlign(interScaledLocalPAM$scaledData,interScaledLocalLPS$scaledData,threshPercent = Thresh)
	
	mdsLPS<-data.frame(mdexp=colMedians(interScaledLocalLPS$scaledData),
										 traj=(1:200)/200,cond="LPS")
	mdsPAM<-data.frame(mdexp=colMedians(interScaledLocalPAM$scaledData),
										 traj=(1:200)/200,cond="PAM")
	
	mdsALL<-c()
	# browser()
	numConn<-length(alignment$align[[1]]$index1)
	for(i in 1:numConn){
		ppam <- alignment$align[[1]]$index1[i]
		plps <- alignment$align[[1]]$index2[i]
		com <- rbind(mdsLPS[plps,],mdsPAM[ppam,])
		com$paired <- i
		mdsALL<-rbind(mdsALL,com)
	}
	
	xpam_aligned<-interScaledLocalPAM$scaledData[,alignment$align[[1]]$index1]
	xlps_aligned<-interScaledLocalLPS$scaledData[,alignment$align[[1]]$index2]
	library(rlist)
	coef_genes<-list()
	for(i in 1:dim(xpam_aligned)[1]){
		coef<-cor.test(xpam_aligned[i,],xlps_aligned[i,],method = "pearson")
		coef_genes<-list.append(coef_genes,coef)
	}
	coef_genes<-unlist(lapply(coef_genes,function(x) return(x$estimate)))
	coef_genes<-unname(unlist(coef_genes))
	df_coef<-data.frame(corr=coef_genes,method=rep(method,length(coef_genes)))
	g<-ggplot(df_coef,aes(x=corr))+geom_density()
	
	
	mdsALL<-mdsALL[c(seq(1,2*numConn,2),seq(2,2*numConn,2)),]
	s<-seq(1,numConn,4)
	mdsALL<-mdsALL[c(s,s+numConn),]
	mdsX <- rbind(mdsLPS,mdsPAM)
	g<-ggplot(mdsALL,aes(x=traj,y=mdexp,group=cond,color=cond))+
		geom_point(size=6)+geom_line(aes(group=paired),size=2,color="grey")+theme_bw(base_size = 40)+
		xlab("Scaled pseudotime")+ylab("Scaled expression")+
		theme(legend.title = element_blank(),
					legend.direction = "horizontal",
					legend.text = element_text(size=60),
					legend.position = c(0.4,0.8))+
		geom_point(data=mdsX,aes(x = traj, y=mdexp, group=cond,color=cond),size=6)
	
	return(list(g,df_coef))
}

compare_correlation<-function(df_coef){
	sample1 <- df_coef[df_coef$method=="Before Integration","corr"]
	methods<-c("DESC", "Harmony", "MNN", "scAlign", "Scanorama", "Seurat V3", "VIPCCA")
	pvalues <- data.frame()
	for(mt in methods){
		print(mt)
		sample2 <- df_coef[df_coef$method==mt,"corr"]
		rho<-t.test(sample1,sample2,paired = T,alternative="two.sided")
		pvalues<- rbind(pvalues, data.frame(pvalue=rho$p.value,method=mt))
	}
	return(pvalues)
}

#glist_alignment<-plot_alignment2("./results/cell_align/desc/")
