Create genome directory first-
To map, follow this

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/newgenome/rawdata/HCT116_siC_1.fq.gz /home/deepa/newgenome/rawdata/HCT116_siC_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/newgenome/star_output/siC/siC1 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/newgenome/rawdata/HCT_siControl_HWN2YCCXX_L5_1.fq.gz /home/deepa/newgenome/rawdata/HCT_siControl_HWN2YCCXX_L5_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/newgenome/star_output/siC/siC2 --outSAMtype BAM SortedByCoordinate 

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/newgenome/rawdata/HCT116_siK_1.fq.gz /home/deepa/newgenome/rawdata/HCT116_siK_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/newgenome/star_output/siT/siT1 --outSAMtype BAM SortedByCoordinate

/home/sb/programfiles/STAR/bin/Linux_x86_64/STAR --runThreadN 18 --genomeDir /home/deepa/newgenome/data/GRCh38/star_indices_overhang100 --readFilesIn /home/deepa/newgenome/rawdata/HCT_siK_HVNYLCCXX_L2_1.fq.gz /home/deepa/newgenome/rawdata/HCT_siK_HVNYLCCXX_L2_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /home/deepa/newgenome/star_output/siT/siT2 --outSAMtype BAM SortedByCoordinate 


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
