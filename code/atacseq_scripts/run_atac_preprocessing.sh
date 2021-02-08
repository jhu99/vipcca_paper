#!/bin/bash
# Filter cell barcodes from BAM file
bam_path="../data/atacseq/10x_pbmc/atac_v1_pbmc_10k_possorted_bam.bam"
cell_path_vipcca="../data/atacseq/10x_pbmc/vipcca.barcode.celltype.txt"
cell_path_seurat="../data/atacseq/10x_pbmc/seurat.barcode.celltype.txt"
cd ./vipcca_bam/
# sinto filterbarcodes -b ${bam_path} -c ${cell_path_vipcca} -p 20
# for f in *.bam;do
# 	SUBSTR=$(echo $f | cut -d'.' -f 1)
# 	echo $SUBSTR
# 	samtools sort ${SUBSTR}.bam > ${SUBSTR}.sort.bam
# 	samtools index ${SUBSTR}.sort.bam
# 	bamCoverage -b ${SUBSTR}.sort.bam -o ${SUBSTR}.bw --binSize 1 --normalizeUsing RPKM --numberOfProcessors 20 --ignoreForNormalization chrX chrY chrM
# done
# 
# multiBigwigSummary BED-file --bwfiles B-cell-progenitor.bw  CD4-Naive.bw Double-negative-T-cell.bw  \
#  pre-B-cell.bw CD14+-Monocytes.bw CD8-effector.bw NK-cell.bw CD16+-Monocytes.bw CD8-Naive.bw \
#  pDC.bw CD4-Memory.bw Dendritic-cell.bw  \
# 	--BED ../data/atacseq/10x_pbmc/markers.bed \
# 	--labels Bpro CD4N DNT Bpre CD14M CD8E NK CD16 CD8N pDC CD4M DC \
# 	-out ./scores_per_transcript.npz --outRawCounts ./scores_per_transcript.tab

multiBigwigSummary bins --bwfiles B-cell-progenitor.bw  CD4-Naive.bw Double-negative-T-cell.bw  \
 pre-B-cell.bw CD14+-Monocytes.bw CD8-effector.bw NK-cell.bw CD16+-Monocytes.bw CD8-Naive.bw \
 pDC.bw CD4-Memory.bw Dendritic-cell.bw  \
	--labels Bpro CD4N DNT Bpre CD14M CD8E NK CD16 CD8N pDC CD4M DC \
	-p 20 -out ./scores_per_transcript_genome.npz --outRawCounts ./scores_per_transcript_genome.tab


cd ../seurat_bam/
# # # sinto filterbarcodes -b ${bam_path} -c ${cell_path_seurat} -p 20
# # for f in *.bam;do
# # 	SUBSTR=$(echo $f | cut -d'.' -f 1)
# # 	echo $SUBSTR
# # 	samtools sort ${SUBSTR}.bam > ${SUBSTR}.sort.bam
# # 	samtools index ${SUBSTR}.sort.bam
# # 	bamCoverage -b ${SUBSTR}.sort.bam -o ${SUBSTR}.bw --binSize 1 --normalizeUsing RPKM --numberOfProcessors 20 --ignoreForNormalization chrX chrY chrM
# # done
# multiBigwigSummary BED-file --bwfiles B-cell-progenitor.bw  CD4-Naive.bw Double-negative-T-cell.bw  \
#  pre-B-cell.bw CD14+-Monocytes.bw CD8-effector.bw NK-cell.bw CD16+-Monocytes.bw CD8-Naive.bw \
#  pDC.bw CD4-Memory.bw Dendritic-cell.bw  \
#  --BED ../data/atacseq/10x_pbmc/markers.bed \
#  --labels Bpro CD4N DNT Bpre CD14M CD8E NK CD16 CD8N pDC CD4M DC \
#  -out ./scores_per_transcript.npz --outRawCounts ./scores_per_transcript.tab

multiBigwigSummary bins --bwfiles B-cell-progenitor.bw  CD4-Naive.bw Double-negative-T-cell.bw  \
 pre-B-cell.bw CD14+-Monocytes.bw CD8-effector.bw NK-cell.bw CD16+-Monocytes.bw CD8-Naive.bw \
 pDC.bw CD4-Memory.bw Dendritic-cell.bw  \
 --labels Bpro CD4N DNT Bpre CD14M CD8E NK CD16 CD8N pDC CD4M DC \
 -p 20 -out ./scores_per_transcript_genome.npz --outRawCounts ./scores_per_transcript_genome.tab

cd ../data/atacseq/bulk_atac_pbmc
for f in *.bam;do
	SUBSTR=$(echo $f | cut -d'.' -f 1)
	echo $SUBSTR
	#samtools sort ${SUBSTR}.bam > ${SUBSTR}.sort.bam
	#samtools index ${SUBSTR}.sort.bam
	#bamCoverage -b ${SUBSTR}.sort.bam -o ${SUBSTR}.bw --binSize 1 --normalizeUsing RPKM --numberOfProcessors 20 --ignoreForNormalization chrX chrY chrM
done
# multiBigwigSummary BED-file --bwfile Bcell.bw CD4.bw CD8.bw NKcell.bw mono.bw \
# 	--BED ../../../data/atacseq/10x_pbmc/markers.bed \
# 	--labels Bcell CD4T CD8T NK MONO \
# 	-out ./scores_per_transcript_marker.npz --outRawCounts ./scores_per_transcript_marker.tab	
# multiBigwigSummary bins --bwfile Bcell.bw CD4.bw CD8.bw NKcell.bw mono.bw \
# 	--labels Bcell CD4T CD8T NK MONO \
# 	-p 20 -out ./scores_per_transcript_genome.npz --outRawCounts ./scores_per_transcript_genome.tab	
# 
# 	
	