computeMatrix reference-point -S ESC.Insulation.10kb.local.bw 2CLC.Insulation.10kb.local.bw -R Up_promoter.bed -b 500000 -a 500000 -o Up_promoter.gz --referencePoint center --binSize 10000 -p 20
plotProfile -m Up_promoter.gz -out Up_promoter.pdf --perGroup --samplesLabel ESC 2CLC --regionsLabel Up_promoter --plotHeight 6 --plotWidth 6 --colors blue red


computeMatrix reference-point -S ESC.Insulation.10kb.local.bw 2CLC.Insulation.10kb.local.bw -R Down_promoter.bed -b 500000 -a 500000 -o Down_promoter.gz --referencePoint center --binSize 10000 -p 20
plotProfile -m Down_promoter.gz -out Down_promoter.pdf --perGroup --samplesLabel ESC 2CLC --regionsLabel Down_promoter --plotHeight 6 --plotWidth 6 --colors blue red


computeMatrix reference-point -S ESC.Insulation.10kb.local.bw 2CLC.Insulation.10kb.local.bw -R MERVL-int-MT2-Mm-expressed.bed -b 500000 -a 500000 -o MERVL-int-MT2-Mm-expressed.gz --referencePoint center --binSize 10000 -p 20
plotProfile -m MERVL-int-MT2-Mm-expressed.gz -out MERVL-int-MT2-Mm-expressed.pdf --perGroup --samplesLabel ESC 2CLC --regionsLabel MERVL-int-MT2-Mm-expressed --plotHeight 6 --plotWidth 6 --colors blue red


computeMatrix reference-point -S ESC.Insulation.10kb.local.bw 2CLC.Insulation.10kb.local.bw -R mm9_promoter_2k2k_strand.bed -b 500000 -a 500000 -o mm9_promoter_2k2k_strand.gz --referencePoint center --binSize 10000 -p 20
plotProfile -m mm9_promoter_2k2k_strand.gz -out mm9_promoter_2k2k_strand.pdf --perGroup --samplesLabel ESC 2CLC --regionsLabel mm9_promoter_2k2k_strand --plotHeight 6 --plotWidth 6 --colors blue red


computeMatrix reference-point -S ESC.Insulation.10kb.local.bw 2CLC.Insulation.10kb.local.bw -R mESC_Dixon2012-raw_TADs.bed -b 500000 -a 500000 -o mESC_Dixon2012-raw_TADs.gz --referencePoint TES --binSize 10000 -p 20
plotProfile -m mESC_Dixon2012-raw_TADs.gz -out mESC_Dixon2012-raw_TADs.pdf --perGroup --samplesLabel ESC 2CLC --regionsLabel mESC_Dixon2012-raw_TADs --plotHeight 6 --plotWidth 6 --colors blue red