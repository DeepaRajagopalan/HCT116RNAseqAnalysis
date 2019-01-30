#istall
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread", version = "3.8")

# load rsubreads
library(Rsubread)
options(scipen=999)

data <- featureCounts(c(
"/home/roberto/deepa/novogene/STAR/HCT116_siC_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT116_siC_DMSO_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT116_siC_JQ1_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT116_siIRF7_DMSO_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT116_siIRF7_JQ1_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT116_siK_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT_DMSO_HWN2YCCXX_L2_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT_JQ1_HWN2YCCXX_L5_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT_siControl_HWN2YCCXX_L5_Aligned.sortedByCoord.out.bam",
"/home/roberto/deepa/novogene/STAR/HCT_siK_HVNYLCCXX_L2_Aligned.sortedByCoord.out.bam"
),
annot.ext="/home/roberto/references/gencode.v26.annotation.gtf",
isGTFAnnotationFile=TRUE,
minMQS=10,
strandSpecific=0,
isPairedEnd=TRUE,
autosort=TRUE,
nthreads=15,
GTF.attrType="gene_name"
)

dat=data[[1]] # extract only the matrix of counts
saveRDS(dat,"novogene_counts.rds")
#
dataset = readRDS("novogene_counts.rds")
