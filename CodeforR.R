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

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library(DESeq2)
