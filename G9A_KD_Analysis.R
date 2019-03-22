#Mapping of Data

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAAARAAPEI-209_1.fq.gz /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAAARAAPEI-209_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/G9A_Analysis/output/shC/C1 --outSAMtype BAM SortedByCoordinate

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAABRAAPEI-210_1.fq.gz /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAABRAAPEI-210_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/G9A_Analysis/output/shC/C2 --outSAMtype BAM SortedByCoordinate

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/G9A_Analysis/FCHNK3WBBXX_L7_HKHUMkjvEAAERABPEI-213_1.fq.gz /home/deepa/G9A_Analysis/FCHNK3WBBXX_L7_HKHUMkjvEAAERABPEI-213_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/G9A_Analysis/output/shB/B1 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAAFRAAPEI-214_1.fq.gz /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAAFRAAPEI-214_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/G9A_Analysis/output/shB/B2 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAACRAAPEI-211_1.fq.gz /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAACRAAPEI-211_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/G9A_Analysis/output/shD/D1 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAADRAAPEI-212_1.fq.gz /home/deepa/G9A_Analysis/FCHNVTNBBXX_L2_HKHUMkjvEAADRAAPEI-212_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/G9A_Analysis/output/shD/D2 --outSAMtype BAM SortedByCoordinate 

library(Rsubread)
options(scipen=999)

data <- featureCounts(c(
"/home/deepa/G9A_Analysis/output/shC/C2Aligned.sortedByCoord.out.bam",
"/home/deepa/G9A_Analysis/output/shC/C1Aligned.sortedByCoord.out.bam",
"/home/deepa/G9A_Analysis/output/shB/B1Aligned.sortedByCoord.out.bam",
"/home/deepa/G9A_Analysis/output/shB/B2Aligned.sortedByCoord.out.bam",
"/home/deepa/G9A_Analysis/output/shD/D1Aligned.sortedByCoord.out.bam",
"/home/deepa/G9A_Analysis/output/shD/D2Aligned.sortedByCoord.out.bam"),
annot.ext="/home/deepa/Annotation/gencode.v27.annotation.gtf",
isGTFAnnotationFile=TRUE,
minMQS=10,
strandSpecific=0,
isPairedEnd=TRUE,
autosort=TRUE,
nthreads=15,
GTF.attrType="gene_name"
)

dat=data[[1]] 
saveRDS(dat,"G9A_Counts.rds")

dataset = readRDS("G9A_Counts.rds")
#DESeq2

library(DESeq2)
designG9A=data.frame(replicate=c("1","2","1","2","1","2"),condition=c("shControl","shControl","shG9A_B","shG9A_B","shG9A_D","shG9A_D"))
dds=DESeqDataSetFromMatrix(countData=dataset, colData=designG9A, design= ~replicate+condition)
#next part is to run the actual differential analysis
ddsgenes=DESeq(dds, test="LRT", full= ~replicate+condition, reduced= ~replicate)
#to extract results in a way that makes sense
dds_results=results(ddsgenes, contrast=c("condition", "shControl", "shG9A_B"))
dds_results2=results(ddsgenes, contrast=c("condition", "shControl", "shG9A_D"))
write.csv(dds_results,"G9A_B_allgenes")
write.csv(dds_results2,"G9A_D_allgenes")

#PlotPCA
ddsgenes_vst <- varianceStabilizingTransformation(ddsgenes)
plotPCA(ddsgenes_vst,ntop=30000,intgroup=c("replicate","condition"))

#heatmap

library(gplots)
library(factoextra)
library(RColorBrewer)
dds_results=results(ddsgenes, contrast=c("condition", "shControl", "shG9A_B"))
dds_results2=results(ddsgenes, contrast=c("condition", "shControl", "shG9A_D"))

dds_genes_heatmap=as.matrix(assay(ddsgenes_vst) [which(abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05 & abs(dds_results2$log2FoldChange)>1 & dds_results2$padj<0.05 ),])
colnames(dds_genes_heatmap)=c("shC1","shC2","shG9A_B1","shG9A_B2","shG9A_D1","shG9A_D2")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(dds_genes_heatmap,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Differentially expressed genes",key.title="Gene expression",cexCol=.8)

interesting_gene =c("EHMT2","MYC","IRF4","RELB")

interesting_genes_heatmap=assay(ddsgenes_vst)[rownames(dds_results) %in% interesting_gene,]
colnames(interesting_genes_heatmap)=c("shC1","shC2","shG9A_B1","shG9A_B2","shG9A_D1","shG9A_D2")
heatmap.2(interesting_genes_heatmap,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
xlab="", ylab="Differentially expressed genes",key.title="Gene expression",cexCol=.8)

#to make Volcano plot
library(gplots)
library(ggplot2)

plot(dds_results$log2FoldChange,-log10(dds_results$padj),xlab="log2FoldChange",
              ylab=expression('-Log'[10]*' p adjusted values'),col=alpha("grey",1),pch=20 )

  abline(v=-1,lty = 2,col="grey")
  abline(v=1,lty = 2,col="grey")
  abline(h=-log10(0.05),lty = 2,col="grey")
  points(dds_results$log2FoldChange[abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05],
       -log10(dds_results$padj)[abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05],
      col=alpha("red",1),pch=20)
      
 points(dds_results$log2FoldChange[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
       -log10(dds_results$padj)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
      col=alpha("blue",1),pch=20)
      
      text(dds_results$log2FoldChange[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
       -log10(dds_results$padj)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
     labels=rownames(dds_results)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")])

  legend("topright", paste("shC",":",length(which(dds_results$log2FoldChange>1 & dds_results$padj<0.05))), bty="n") 
  legend("topleft", paste("shG9A_B",":",length(which(dds_results$log2FoldChange<(-1) & dds_results$padj<0.05))), bty="n")
  
points(dds_results$log2FoldChange[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))],
       -log10(dds_results$padj)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))],
      col=alpha("black",1),pch=20)

text(dds_results$log2FoldChange[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))]+0.2, #to move the position of the label by 0.2 spaces.
       -log10(dds_results$padj)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))]+1,
     labels=rownames(dds_results)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))])

#for shG9A_D

library(gplots)
library(ggplot2)

plot(dds_results2$log2FoldChange,-log10(dds_results2$padj),xlab="log2FoldChange",
              ylab=expression('-Log'[10]*' p adjusted values'),col=alpha("grey",1),pch=20 )

  abline(v=-1,lty = 2,col="grey")
  abline(v=1,lty = 2,col="grey")
  abline(h=-log10(0.05),lty = 2,col="grey")
  points(dds_results2$log2FoldChange[abs(dds_results2$log2FoldChange)>1 & dds_results2$padj<0.05],
       -log10(dds_results2$padj)[abs(dds_results2$log2FoldChange)>1 & dds_results2$padj<0.05],
      col=alpha("red",1),pch=20)
      
 points(dds_results2$log2FoldChange[which(rownames(dds_results2)=="EHMT2" | rownames(dds_results2)=="MYC"| rownames(dds_results2)=="IRF4"| rownames(dds_results2)=="RELB")],
       -log10(dds_results2$padj)[which(rownames(dds_results2)=="EHMT2" | rownames(dds_results2)=="MYC"| rownames(dds_results2)=="IRF4"| rownames(dds_results2)=="RELB")],
      col=alpha("blue",1),pch=20)
      
      text(dds_results2$log2FoldChange[which(rownames(dds_results2)=="EHMT2" | rownames(dds_results2)=="MYC"| rownames(dds_results2)=="IRF4"| rownames(dds_results2)=="RELB")],
       -log10(dds_results2$padj)[which(rownames(dds_results2)=="EHMT2" | rownames(dds_results2)=="MYC"| rownames(dds_results2)=="IRF4"| rownames(dds_results2)=="RELB")],
     labels=rownames(dds_results2)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results2)=="MYC"| rownames(dds_results2)=="IRF4"| rownames(dds_results2)=="RELB")])

  legend("topright", paste("shC",":",length(which(dds_results2$log2FoldChange>1 & dds_results2$padj<0.05))), bty="n") 
  legend("topleft", paste("shG9A_D",":",length(which(dds_results2$log2FoldChange<(-1) & dds_results2$padj<0.05))), bty="n")
  
  points(dds_results2$log2FoldChange[c(which.max(dds_results2$log2FoldChange),which.min(dds_results2$log2FoldChange))],
       -log10(dds_results$padj)[c(which.max(dds_results2$log2FoldChange),which.min(dds_results2$log2FoldChange))],
      col=alpha("black",1),pch=20)
      
      text(dds_results2$log2FoldChange[c(which.max(dds_results2$log2FoldChange),which.min(dds_results2$log2FoldChange))]+0.2, #to move the position of the label by 0.2 spaces.
       -log10(dds_results2$padj)[c(which.max(dds_results2$log2FoldChange),which.min(dds_results2$log2FoldChange))]+1,
     labels=rownames(dds_results2)[c(which.max(dds_results2$log2FoldChange),which.min(dds_results2$log2FoldChange))])

#to consider the shRNAs as replicates

library(DESeq2)
dataset = readRDS("G9A_Counts.rds")
designG9A=data.frame(condition=c("shControl","shControl","shG9A","shG9A","shG9A","shG9A"))
dds=DESeqDataSetFromMatrix(countData=dataset, colData=designG9A, design= ~condition)
#next part is to run the actual differential analysis
ddsgenes=DESeq(dds, test="Wald")
#to extract results in a way that makes sense
dds_results=results(ddsgenes, contrast=c("condition", "shControl", "shG9A"))

ddsgenes_vst <- varianceStabilizingTransformation(ddsgenes)
plotPCA(ddsgenes_vst,ntop=30000,intgroup=c("condition"))

#volcano plot

library(gplots)
library(ggplot2)

plot(dds_results$log2FoldChange,-log10(dds_results$padj),xlab="log2FoldChange",
              ylab=expression('-Log'[10]*' p adjusted values'),col=alpha("grey",1),pch=20 )

  abline(v=-1,lty = 2,col="grey")
  abline(v=1,lty = 2,col="grey")
  abline(h=-log10(0.05),lty = 2,col="grey")
  points(dds_results$log2FoldChange[abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05],
       -log10(dds_results$padj)[abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05],
      col=alpha("red",1),pch=20)
      
 points(dds_results$log2FoldChange[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
       -log10(dds_results$padj)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
      col=alpha("blue",1),pch=20)
      
      text(dds_results$log2FoldChange[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
       -log10(dds_results$padj)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")],
     labels=rownames(dds_results)[which(rownames(dds_results)=="EHMT2" | rownames(dds_results)=="MYC"| rownames(dds_results)=="IRF4"| rownames(dds_results)=="RELB")])

  legend("topright", paste("shC",":",length(which(dds_results$log2FoldChange>1 & dds_results$padj<0.05))), bty="n") 
  legend("topleft", paste("shG9A",":",length(which(dds_results$log2FoldChange<(-1) & dds_results$padj<0.05))), bty="n")
  
points(dds_results$log2FoldChange[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))],
       -log10(dds_results$padj)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))],
      col=alpha("black",1),pch=20)

text(dds_results$log2FoldChange[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))]+0.2, #to move the position of the label by 0.2 spaces.
       -log10(dds_results$padj)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))]+1,
     labels=rownames(dds_results)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))])

#heatmap

library(gplots)
library(factoextra)
library(RColorBrewer)

dds_genes_heatmap=as.matrix(assay(ddsgenes_vst) [which(abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05),])
colnames(dds_genes_heatmap)=c("shC1","shC2","shG9A","shG9A","shG9A","shG9A")
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
heatmap.2(dds_genes_heatmap,col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="Differentially expressed genes",key.title="Gene expression",cexCol=.8)

rownames(dds_results)[dds_results$log2FoldChange>1 & dds_results$padj<0.05]
write.csv(rownames(dds_results)[dds_results$log2FoldChange>1 & dds_results$padj<0.05],"downregulated_genes_shG9A")
