# for Hi-C data analysis, raw reads were trimmed to 50 bp and perform analysis using HiC-Pro
# process sisHi-C data
HiC-Pro -i rawdata/ -o HICPRO1 -c config-hicpro.txt

# Then further normalize matrix data using HiC-Pro-normalize.R for different resolution
Rscript HiC-Pro-normalize.R ESC_500000_iced.matrix
Rscript HiC-Pro-normalize.R 2CLC_500000_iced.matrix

# convert HiC-Pro data to .hic format to visualize in Juicer (Fig 1A)
awk '{$4=$4!="+"; $7=$7!="+"} $2<=$5{print $1, $4, $2, $3, 0, $7, $5, $6, 1, $11, $12 }$5<$2{ print $1, $7, $5, $6, 0, $4, $2, $3, 1, $12, $11 }' ESC.allValidPairs | sort  -k3,3d  -k7,7d -S 90 > ESC_allValidPairs.pre_juicebox_sorted
awk '{$4=$4!="+"; $7=$7!="+"} $2<=$5{print $1, $4, $2, $3, 0, $7, $5, $6, 1, $11, $12 }$5<$2{ print $1, $7, $5, $6, 0, $4, $2, $3, 1, $12, $11 }' 2CLC.allValidPairs | sort  -k3,3d  -k7,7d -S 90 > 2CLC_allValidPairs.pre_juicebox_sorted
java -jar juicer_tools_1.22.01.jar pre -r 5000,10000,20000,25000,50000,100000,250000,500000 -k VC,VC_SQRT,KR -d 2CLC_allValidPairs.pre_juicebox_sorted 2CLC_allValidPairs.hic mm9
java -jar juicer_tools_1.22.01.jar pre -r 5000,10000,20000,25000,50000,100000,250000,500000 -k VC,VC_SQRT,KR -d ESC_allValidPairs.pre_juicebox_sorted ESC_allValidPairs.hic mm9

# convert 500kb resolution HiC-Pro matrix data (sparse) to Dense data using sparse2dense.py for each chromosome
python sparse2dense.py ESC_500000_abs.bed ESC_500000_iced.matrix.norm.matrix
python sparse2dense.py 2CLC_500000_abs.bed 2CLC_500000_iced.matrix.norm.matrix

# calculate PC1 values using cworld (Dekker lab)
for file in *dense.matrix.txt
do
perl /nethome/yezhang_zhu/tools/cworld-dekker-master/scripts/perl/matrix2compartment.pl -i $file 
done

# combine all PC1 bedGraph files for each chromosome together and then visualize in IGV (Fig 1B, Fig S2E)
cat *bedGraph | grep -v "track" > ESC.PC1.bedGraph
cat *bedGraph | grep -v "track" > 2CLC.PC1.bedGraph

sed 's/nan/0/g' ESC.PC1.bedGraph > ESC.PC1.nanto0.bedGraph
sed 's/nan/0/g' 2CLC.PC1.bedGraph > 2CLC.PC1.nanto0.bedGraph

# reorder the bins according to their PC1 values using reorder-bins-according-to-PC1-values.R

# plot PC1 values using plotPC1.R (Fig S2A upper)

# calculate compartment heatmap using Compartmentheatmap.py and then plot using plotHeatmap.R (Fig S2A lower)
python Compartmentheatmap.py ESC_500000_abs.bed ESC_500000_iced.matrix.norm.matrix Order.Compartment.ESC.bed ESC.compartment.txt
python Compartmentheatmap.py 2CLC_500000_abs.bed 2CLC_500000_iced.matrix.norm.matrix Order.Compartment.2CLC.bed 2CLC.compartment.txt

# calculate compartment strength for each chromosome and average strength using Compartmentstrength.py and then plot using Fig1C.R (Fig 1C), FigS2B.R (Fig S2B)
python Compartmentstrength.py ESC_500000_abs.bed 4w-2Cn_500000_iced.matrix.norm.matrix Order.Compartment.4w-2Cn.bed 4w-2Cn.compartment.5.txt > 4w-2Cn.compartment.strength.txt
python Compartmentstrength.py ESC_500000_abs.bed 4w-2Cp_500000_iced.matrix.norm.matrix Order.Compartment.4w-2Cp.bed 4w-2Cp.compartment.5.txt > 4w-2Cp.compartment.strength.txt
python Compartmentstrength.py ESC_500000_abs.bed 10wESC_500000_iced.matrix.norm.matrix Order.Compartment.10wESC.bed 10wESC.compartment.5.txt > 10wESC.compartment.strength.txt
python Compartmentstrength.py ESC_500000_abs.bed 10w2CLC_500000_iced.matrix.norm.matrix Order.Compartment.10w2CLC.bed 10w2CLC.compartment.5.txt > 10w2CLC.compartment.strength.txt

# calculate genmic intervals that undergo A to B compartment of B to A compartment changes using FigS2C.R (Fig S2C)

# obtain B to A compartment covered genes and then overlap with Up-regulated genes detected in RNA-seq using FigS2DFG.R
bedtools intersect -u -wa -a mm9_TSS_GeneName.bed -b BtoA.bed > mm9_TSS_GeneName_BtoA.bed

# calculate enhancer enrichment in A to B or B to A bins using overlapOSN-with-BtoA-AtoB.sh
rm enhancer.enrichment.txt
for i in {1..100}
do
bash overlapOSN-with-BtoA-AtoB.sh >> enhancer.enrichment.txt
done 
# and then plot using FigS2DFG.R

# calculate expression level changes in RNA-seq using calculate-RNA-expression-changes.sh and then plot using FigS2DFG.R

# convert 20kb resolution HiC-Pro matrix data (sparse) to Dense data using sparse2dense.py for each chromosome
python sparse2dense.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix
python sparse2dense.py 2CLC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix

# calculate insulation scores and then call TADs using cworld (Dekker lab)
for file in *dense.matrix.txt
do
perl /nethome/yezhang_zhu/tools/cworld-dekker-master/scripts/perl/matrix2insulation.pl -i $file 
done

for file in *insulation; do perl /nethome/yezhang_zhu/tools/cworld-dekker-master/scripts/perl/insulation2tads.pl -i $file -b ${file}.boundaries; done

# combine all TAD files for each chromosome together and then calculate overlap between ESC and 2CLC and make venn Diagram (Fig S2I)
cat *tads.bed | grep -v "track" > ESC.tad.bed
cat *tads.bed | grep -v "track" > 2CLC.tad.bed
awk '{if($3-$2 >= 100000) print $0}' ESC.tad.bed > ESC.tad.bed.2
awk '{if($3-$2 >= 100000) print $0}' 2CLC.tad.bed > 2CLC.tad.bed.2
mv ESC.tad.bed.2 ESC.tad.bed
mv 2CLC.tad.bed.2 2CLC.tad.bed
bedtools intersect -wa -wb -f 0.75 -r -a ESC.tad.bed -b 2CLC.tad.bed > stable.tad.bed
bedtools intersect -v -wa -wb -f 0.75 -r -a ESC.tad.bed -b stable.tad.bed > ESC.nostable.bed
bedtools intersect -v -wb -wb -f 0.75 -r -a 2CLC.tad.bed -b stable.tad.bed > 2CLC.nostable.bed


# convert 10kb resolution HiC-Pro matrix data (sparse) to Dense data using sparse2dense.py for each chromosome
python sparse2dense.py ESC_10000_abs.bed ESC_10000_iced.matrix.norm.matrix
python sparse2dense.py 2CLC_10000_abs.bed 2CLC_10000_iced.matrix.norm.matrix

# calculate insulation scores and then normalize using Crane-local.R (normalize using local neighboring bins)
for file in *dense.matrix.txt
do
perl /nethome/yezhang_zhu/tools/cworld-dekker-master/scripts/perl/matrix2insulation.pl -i $file 
done
cd 2CLC
Rscript Crane-local.R
cat *out > 2CLC.Insulation.10kb.bedGraph

cd ESC
Rscript Crane-local.R
cat *out > ESC.Insulation.10kb.bedGraph

grep -v "NA" ESC.Insulation.10kb.bedGraph > ESC.Insulation.10kb.local.bedGraph
grep -v "NA" 2CLC.Insulation.10kb.bedGraph > 2CLC.Insulation.10kb.local.bedGraph

awk -v OFS="\t" '{print $1,$2,$3,-1*$4}' ESC.Insulation.10kb.local.bedGraph > ESC.Insulation.10kb.local.local.bedGraph
awk -v OFS="\t" '{print $1,$2,$3,-1*$4}' 2CLC.Insulation.10kb.local.bedGraph > 2CLC.Insulation.10kb.local.local.bedGraph

wigToBigWig -clip ESC.Insulation.10kb.local.local.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes ESC.Insulation.10kb.local.bw
wigToBigWig -clip 2CLC.Insulation.10kb.local.local.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes 2CLC.Insulation.10kb.local.bw

# generate average profile figures using generate-insulation-profile.sh (Fig S2 H, J, K)

gunzip mESC_Dixon2012-raw_TADs.gz
# further compare insulation scores at TAD boundaries using Fig1DE.R (Fig 1D, E)

# extract OE matrix at specific genomic loci in Fig 1F and Fig 3H using Extract.OE.matrix.py
python Extract.OE.matrix.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix 2CLC_20000_iced.matrix.norm.matrix chr16:31000000-32400000 ESC.Fig1F 2CLC.Fig1F
python Extract.OE.matrix.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix 2CLC_20000_iced.matrix.norm.matrix chr5:33200000-34800000 ESC.Fig3H 2CLC.Fig3H
# and then plot using Fig1F3H using Fig1F3H.R

# calculate Directionality index using DI.py
python DI.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix
python DI.py ESC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix
python DI.py ESC_20000_abs.bed Untreat_20000_iced.matrix.norm.matrix
python DI.py ESC_20000_abs.bed auxin-2days_20000_iced.matrix.norm.matrix
python DI.py ESC_20000_abs.bed washoff-2days_20000_iced.matrix.norm.matrix

cat *bedGraph > 2CLC.DI.bedGraph
sed 's/NA/0/g' 2CLC.DI.bedGraph > 2CLC.DI.NAto0.bedGraph
wigToBigWig -clip 2CLC.DI.NAto0.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes 2CLC.DI.NAto0.bw

cat *bedGraph > ESC.DI.bedGraph
sed 's/NA/0/g' ESC.DI.bedGraph > ESC.DI.NAto0.bedGraph
wigToBigWig -clip ESC.DI.NAto0.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes ESC.DI.NAto0.bw

cat *bedGraph > Untreat.DI.bedGraph
sed 's/NA/0/g' Untreat.DI.bedGraph > Untreat.DI.NAto0.bedGraph
wigToBigWig -clip Untreat.DI.NAto0.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes Untreat.DI.NAto0.bw

cat *bedGraph > auxin-2days.DI.bedGraph
sed 's/NA/0/g' auxin-2days.DI.bedGraph > auxin-2days.DI.NAto0.bedGraph
wigToBigWig -clip auxin-2days.DI.NAto0.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes auxin-2days.DI.NAto0.bw

cat *bedGraph > washoff-2days.DI.bedGraph
sed 's/NA/0/g' washoff-2days.DI.bedGraph > washoff-2days.DI.NAto0.bedGraph
wigToBigWig -clip washoff-2days.DI.NAto0.bedGraph ~/tools/ucsc_tools/mm9.chrom.sizes washoff-2days.DI.NAto0.bw

# then generate average plot figures using Fig1G4E.sh and visualize Directionality index in IGV (Fig 1F lower)

# extract raw matrix at specific genomic loci in Fig S3A using Extract.raw.matrix.py
python Extract.raw.matrix.py ESC_10000_abs.bed ESC_10000_iced.matrix.norm.matrix 2CLC_10000_iced.matrix.norm.matrix chr6:122630000-122740000 ESC.Nanog 2CLC.Nanog
python Extract.raw.matrix.py ESC_10000_abs.bed ESC_10000_iced.matrix.norm.matrix 2CLC_10000_iced.matrix.norm.matrix chr3:34470000-34740000 ESC.Sox2 2CLC.Sox2
python Extract.raw.matrix.py ESC_10000_abs.bed ESC_10000_iced.matrix.norm.matrix 2CLC_10000_iced.matrix.norm.matrix chr4:55460000-55570000 ESC.Klf4 2CLC.Klf4
# and then plot using Fig S3A using FigS3A.R and extract the interaction score of indicated EPI and plot using FigS3B.R

# calculate PS curve using PScurve1.py and PScurve2.py
python PScurve1.py ESC_5000_abs.bed 4w-2Cn-facs_5000_iced.matrix.norm.matrix > 4w-2Cn-facs.PS.txt
python PScurve1.py ESC_5000_abs.bed 4w-2Cp-facs_5000_iced.matrix.norm.matrix > 4w-2Cp-facs.PS.txt
python PScurve1.py ESC_5000_abs.bed 10wESC-sisHiC-r2-yjl_5000_iced.matrix.norm.matrix > 10wESC.PS.txt
python PScurve1.py ESC_5000_abs.bed 10w2CLC-sisHiC-r2-yjl_5000_iced.matrix.norm.matrix > 10w2CLC.PS.txt

python PScurve2.py 4w-2Cn-facs.PS.txt > 4w-2Cn-facs.PS.2.txt
python PScurve2.py 4w-2Cp-facs.PS.txt > 4w-2Cp-facs.PS.2.txt
python PScurve2.py 10wESC.PS.txt > 10wESC.PS.2.txt
python PScurve2.py 10w2CLC.PS.txt > 10w2CLC.PS.2.txt
# and then plot using Fig1H.R (Fig 1H)

# calculate heatmap around TAD boundaries using Boundary.flank.py and plot using Boundary.flank.R (Fig 1I)
python Boundary.flank.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix mESC_Dixon2012-raw_TADs.bed ESC.boundary.txt 
python Boundary.flank.py ESC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix mESC_Dixon2012-raw_TADs.bed 2CLC.boundary.txt 

# calculate TAD strength using TADstrength.py and plot using TADstrength.R (Fig 1J) [also used in Fig S2L]
python TADstrength.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix mESC_Dixon2012-raw_TADs.bed 20000 > ESC.tadstrength.txt 
python TADstrength.py ESC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix mESC_Dixon2012-raw_TADs.bed 20000 > 2CLC.tadstrength.txt 

# calculate inter TAD frequency using InterTAD.py and plot using InterTAD.R (Fig 1K)
python InterTAD.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix mESC_Dixon2012-raw_TADs.bed 20000 > ESC.interTAD.txt
python InterTAD.py ESC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix mESC_Dixon2012-raw_TADs.bed 20000 > 2CLC.interTAD.txt

# calculate loop-loop heatmap using Loop.aggregate.py and plot using Loop.aggregate.R (Fig 2A, 4F, S2M)
python Loop.aggregate.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt ESC.Loop.txt 
python Loop.aggregate.py ESC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt 2CLC.Loop.txt 
python Loop.aggregate.py ESC_20000_abs.bed Untreat_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt Untreat.Loop.txt 
python Loop.aggregate.py ESC_20000_abs.bed auxin-2days_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt auxin-2days.Loop.txt 
python Loop.aggregate.py ESC_20000_abs.bed washoff-2days_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt washoff-2days.Loop.txt 
python Loop.aggregate.py ESC_20000_abs.bed early-2Cell_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt early-2Cell.Loop.txt 
python Loop.aggregate.py ESC_20000_abs.bed ICM_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt ICM.Loop.txt 

# calculate loop-loop strength using loopstrength.py and plot using loopstrength.R (Fig 2B, S2M) 
python loopstrength.py ESC_20000_abs.bed ESC_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt > ESC.Loopstrength.txt 
python loopstrength.py ESC_20000_abs.bed 2CLC_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt > 2CLC.Loopstrength.txt 
python loopstrength.py ESC_20000_abs.bed early-2Cell_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt > early2Cell.Loopstrength.txt 
python loopstrength.py ESC_20000_abs.bed ICM_20000_iced.matrix.norm.matrix GSE63525_CH12-LX_HiCCUPS_looplist.txt > ICM.Loopstrength.txt 


# calculate EPI heatmap using EPI.aggregate.py and plot using EPI.aggregate.R (Fig 2C, 4G, S2N)
python EPI.aggregate.py ESC_10000_abs.bed ESC_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed ESC.EPI.txt 
python EPI.aggregate.py ESC_10000_abs.bed 2CLC_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed 2CLC.EPI.txt 
python EPI.aggregate.py ESC_10000_abs.bed early-2Cell_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed early-2Cell.EPI.txt 
python EPI.aggregate.py ESC_10000_abs.bed ICM_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed ICM.EPI.txt 
python EPI.aggregate.py ESC_10000_abs.bed Untreat_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed Untreat.EPI.txt 
python EPI.aggregate.py ESC_10000_abs.bed auxin-2days_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed auxin-2days.EPI.txt 
python EPI.aggregate.py ESC_10000_abs.bed washoff-2days_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed washoff-2days.EPI.txt 


# calculate EPI strength using EPIstrength.py and plot using EPIstrength.R (Fig 2D, S2N) 
python EPIstrength.py ESC_10000_abs.bed ESC_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed > ESC.EPIstrength.txt 
python EPIstrength.py ESC_10000_abs.bed 2CLC_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed > 2CLC.EPIstrength.txt 
python EPIstrength.py ESC_10000_abs.bed early-2Cell_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed > early2Cell.EPIstrength.txt 
python EPIstrength.py ESC_10000_abs.bed ICM_10000_iced.matrix.norm.matrix RUAN_enhancer_promoter.bed > ICM.EPIstrength.txt 
