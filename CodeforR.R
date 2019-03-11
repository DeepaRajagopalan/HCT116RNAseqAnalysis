Create genome directory first-
To map, follow this

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/Hep3B_Analysis/DMSO_1.fq.gz /home/deepa/Hep3B_Analysis/DMSO_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/Hep3B_Analysis/star_output/DMSO/D1 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/Hep3B_Analysis/Hep3BDMSO_HWN2YCCXX_L5_1.clean.fq.gz  /home/deepa/Hep3B_Analysis/Hep3BDMSO_HWN2YCCXX_L5_2.clean.fq.gz  --readFilesCommand zcat --outFileNamePrefix /home/deepa/Hep3B_Analysis/star_output/DMSO/D2 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/Hep3B_Analysis/JQ1_1.fq.gz /home/deepa/Hep3B_Analysis/JQ1_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/Hep3B_Analysis/star_output/JQ1/J1 --outSAMtype BAM SortedByCoordinate

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/Hep3B_Analysis/Hep3BJQ1_HWN2YCCXX_L5_1.clean.fq.gz  /home/deepa/Hep3B_Analysis/Hep3BJQ1_HWN2YCCXX_L5_2.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/Hep3B_Analysis/star_output/JQ1/J2 --outSAMtype BAM SortedByCoordinate 


# NOT Necessary to use
#samtools view -bS Aligned.out.sam > Aligned.out.bam
#samtools sort Aligned.out.bam > sorted.bam
#samtools index -b sorted.bam


#istall for R version 3.5
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread", version = "3.8")

#to install Rsubreads  for R3.3.3.1 use this
source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")

# load rsubreads
library(Rsubread)
options(scipen=999)

data <- featureCounts(c(
"/home/deepa/newgenome/star_output/siC/siC2Aligned.sortedByCoord.out.bam",
"/home/deepa/newgenome/star_output/siC/siC1Aligned.sortedByCoord.out.bam",
"/home/deepa/newgenome/star_output/siT/siT1Aligned.sortedByCoord.out.bam",
"/home/deepa/newgenome/star_output/siT/siT2Aligned.sortedByCoord.out.bam"),
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
saveRDS(dat,"HCT116_counts.rds")
#
dataset = readRDS("HCT116_counts.rds")
#
dat=data[[1]]  # extract only the matrix of counts

# for the DEseq2 analysis
# first install DESeq2 package)
#Don't run on computer R if DESeq is already installed
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library(DESeq2)
designHCT116=data.frame(replicate=c("2","1","1","2"),condition=c("siControl","siControl","siTIP60","siTIP60"))
dds=DESeqDataSetFromMatrix(countData=dataset, colData=designHCT116, design= ~replicate+condition)
#next part is to run the actual differential analysis
ddsgenes=DESeq(dds, test="LRT", full= ~replicate+condition, reduced= ~replicate)
#to extract results in a way that makes sense
dds_results=results(ddsgenes, contrast=c("condition", "siControl", "siTIP60"))
write.csv(dds_results,"HCT_siC_vs_siK_allgenes")

#to count the number of significantly expressed genes

#for PCA
ddsgenes_vst <- varianceStabilizingTransformation(ddsgenes) #normalize by library size and transform to something like log 2 (taking into consideration base mean, penalizes genes that are lowly expressed)

plotPCA(ddsgenes_vst,ntop=30000,intgroup=c("replicate","condition"))
#tomake Volcano plot

plot(dds_results$log2FoldChange,-log10(dds_results$padj),xlab="log2FoldChange",
              ylab=expression('-Log'[10]*' p adjusted values'),col=alpha("grey",1),pch=20 )

  abline(v=-1,lty = 2,col="grey")
  abline(v=1,lty = 2,col="grey")
  abline(h=-log10(0.05),lty = 2,col="grey")
  points(dds_results$log2FoldChange[abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05],
       -log10(dds_results$padj)[abs(dds_results$log2FoldChange)>1 & dds_results$padj<0.05],
      col=alpha("red",1),pch=20)

  legend("topright", paste("siC",":",length(which(dds_results$log2FoldChange>1 & dds_results$padj<0.05))), bty="n") 
  legend("topleft", paste("siK",":",length(which(dds_results$log2FoldChange<(-1) & dds_results$padj<0.05))), bty="n") 

#to highlight specific genes
points(dds_results$log2FoldChange[which(rownames(dds_results)=="KAT5" | rownames(dds_results)=="IRF7")],
       -log10(dds_results$padj)[which(rownames(dds_results)=="KAT5" | rownames(dds_results)=="IRF7")],
      col=alpha("blue",1),pch=20)
#to extract the location of desired genes from the rows
which(rownames(dds_results)=="KAT5" | rownames(dds_results)=="IRF7")

text(dds_results$log2FoldChange[which(rownames(dds_results)=="KAT5" | rownames(dds_results)=="IRF7")],
       -log10(dds_results$padj)[which(rownames(dds_results)=="KAT5" | rownames(dds_results)=="IRF7")],
     labels=rownames(dds_results)[which(rownames(dds_results)=="KAT5" | rownames(dds_results)=="IRF7")])
#to extract the most different genes  and then plot
which.max(dds_results$log2FoldChange) | which.min(dds_results$log2FoldChange)
points(dds_results$log2FoldChange[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))],
       -log10(dds_results$padj)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))],
      col=alpha("black",1),pch=20)

text(dds_results$log2FoldChange[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))]+0.2, #to move the position of the label by 0.2 spaces.
       -log10(dds_results$padj)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))]+1,
     labels=rownames(dds_results)[c(which.max(dds_results$log2FoldChange),which.min(dds_results$log2FoldChange))])



