mkdir bw
mkdir bedgraph
for file in Neg.bam Pos.bam
do
 readNum=` samtools view ${file} | wc -l `
 scaleNum=`echo "scale=6;1000000 / $readNum;" | bc`
 bedtools genomecov -bg -split -scale $scaleNum -ibam ${file} -g ~/tools/ucsc_tools/mm9.chrom.sizes > ./bedgraph/${file/.bam/}.bdg
 bedGraphToBigWig ./bedgraph/${file/.bam/}.bdg ~/tools/ucsc_tools/mm9.chrom.sizes ./bw/${file/.bam/}.bw
done

