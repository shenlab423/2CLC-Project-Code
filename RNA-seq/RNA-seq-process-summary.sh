# process ESC and 2CLC RNA-seq data

GTF=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf
GENOME_FA=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome.fa
GENOME=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

for file in `ls *R1*.fastq`
do 
{
trim-c $file ${file%%_*}_1.fastq 0 50
} &
done
wait

for file in `ls *_1.fastq`
do
outdir=${file%_*}_tophat
tophat -o ${outdir} -p 20 -G ${GTF} ${GENOME} ${file} 
done

for filedir in `ls -F | grep 'tophat/$'`
do
{
bash filter_unique_reads_from_tophat.sh ${filedir}accepted_hits.bam ${filedir}accepted_hits_unique
}&
done
wait

# convert bam files to bw using visualize.sh and then visualize in IGV (Fig S2E)
# obtain count matrix using htseq-count
for file in *bam; do htseq-count -s no -t exon -f bam $file mm9.genes.exon.fix.gtf > ${file/.bam/}.count.txt & done

# combine count data and then analyzed using edgeR package (Fig2H.R)
paste *count.txt | cut -f 1,2,4,6,8,10,12,14,16 > Shenlab.2C.exon.fix.2.count

# identify DEGs using FigS1EH.R

# exam gene expression of DEGs in mouse embryos using FigS1EH.R

# perform pathway analysis using IPA and plot using FigS1IJ.R

# analysis super-enhancer neighboring genes expression using Fig2G.R

# generate GSEA figures using gsea-3.0 with SuperEnhancer-Nearest-gene_sets.gmt GSEA gene set

# calculate expression level of repeats using repeat_expression.sh
paste *ESC*rpkm > Neg.RNA.repeats.txt
paste *2CLC*rpkm > Pos.RNA.repeats.txt

# and then analyze using FigS1FG.R

# process ESC and 2CLC RNA-seq data in CTCF-AID cells

GTF=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf
GENOME_FA=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome.fa
GENOME=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

for file in `ls *R1*.fastq`
do 
{
trim-c $file ${file%%_*}_1.fastq 0 50
} &
done
wait

for file in `ls *_1.fastq`
do
outdir=${file%_*}_tophat
tophat -o ${outdir} -p 20 -G ${GTF} ${GENOME} ${file} #--no-novel-juncs 
done

for filedir in `ls -F | grep 'tophat/$'`
do
{
bash filter_unique_reads_from_tophat.sh ${filedir}accepted_hits.bam ${filedir}accepted_hits_unique
}&
done
wait

for file in *bam; do htseq-count -s no -t exon -f bam $file mm9.genes.exon.fix.gtf > ${file/.bam/}.count.txt & done
# combine count data and then analyzed using edgeR package (CTCF-AID-expression.R)
paste *count.txt > RNA-IAA-2d-2C.count.txt

# compare ESC and 2CLC RNA-seq data in WT and CTCF-AID cells using FigS5E.R

# compare architecture gene expression in ESC and 2CLC RNA-seq data of WT and CTCF-AID cells using FigS5F.R

# processing published CTCF-AUD treat and RAD21-AID treat RNA-seq data using CTCF-AID-RAD21-AID-RNA-seq-process.sh

# compare the fold change of 2CLC-specific genes, all genes, and pluripotent genes in CTCF-AID and RAD21-AID cells using Fig4CS5C.R

