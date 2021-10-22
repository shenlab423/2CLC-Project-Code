# process ATAC-seq, ChIP-seq for histone modifications in ESC, 2CLC and mouse embryos using ATAC-histone-modification-process.sh
# and then visualize in IGV (Fig 3H, S3A)

# calculate ESC and 2CLC-specific peaks
bash cal_specific_peaks.sh ESC-H3K27ac_peaks.broadPeak 2CLC-H3K27ac_peaks.broadPeak H3K27ac ESC-H3K27ac.sam 2CLC-H3K27ac.sam
bash cal_specific_peaks.sh ESC-H3K27me3_peaks.broadPeak 2CLC-H3K27me3_peaks.broadPeak H3K27me3 ESC-H3K27me3.sam 2CLC-H3K27me3.sam
bash cal_specific_peaks.sh ESC-H3K4me3_peaks.broadPeak 2CLC-H3K4me3_peaks.broadPeak H3K4me3 ESC-H3K4me3.sam 2CLC-H3K4me3.sam
bash cal_specific_peaks.sh ESC-H3K4me1_peaks.broadPeak 2CLC-H3K4me1_peaks.broadPeak H3K4me1 ESC-H3K4me1.sam 2CLC-H3K4me1.sam
bash cal_specific_peaks.sh ESC-H3K9me3_peaks.broadPeak 2CLC-H3K9me3_peaks.broadPeak H3K9me3 ESC-H3K9me3.sam 2CLC-H3K9me3.sam
# and then make venn Diagram in Fig S3C using overlap.R

# test the distribution of peaks using ceasBW and collect all data into Distribution.csv
ceasBW -w ESC-H3K27ac.bw -b H3K27ac_ESC.bed -g ~/tools/mm9.refGene -l ~/tools/mm9.chrom.sizes
# plot the distribution using Distribution.R (Fig S3D)

# annotate H3K27ac peaks in 2CLC-specific
cut -f 1,2,3 H3K27ac_2C_UP_rpkm.bdg > H3K27ac_2C_UP_rpkm.bed
awk -v OFS="\t" '{if($5>4*$4) print $1,$2,$3}' H3K27ac_2C_UP_rpkm.bdg > H3K27ac_2C_UP_rpkm_FC4.bed
annotatePeaks.pl H3K27ac_2C_UP_rpkm.bed none -gtf /nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf > H3K27ac_2C_UP_rpkm.bed.annotated.txt
annotatePeaks.pl H3K27ac_2C_UP_rpkm_FC4.bed none -gtf /nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf > H3K27ac_2C_UP_rpkm_FC4.bed.annotated.txt

# obtain all putative 2C enhancers and putative 2C enhancers (fold change > 4) using Get-putative-enhancers.R

# generate Fig 3A,B,C,D,S4A using deeptools(computeMatrix, and plotProfile)

# process DUX ChIP-seq data using DUX-ChIP-process.sh and then generate Fig 3F using deeptools(computeMatrix, and plotProfile)

# process CTCF CUT&TAG data using CTCF-CUTTAG-process.sh and then generate Fig S5G using deeptools(computeMatrix, and plotProfile)