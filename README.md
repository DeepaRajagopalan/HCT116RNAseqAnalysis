#HCT116RNAseqAnalysis

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

