computeMatrix reference-point -S ESC.DI.NAto0.bw 2CLC.DI.NAto0.bw -R mESC_Dixon2012-raw_TADs.boundary.bed -b 600000 -a 600000 -o Boundary.gz --referencePoint center --binSize 20000 -p 20

plotProfile -m Boundary.gz -out Boundary.pdf --perGroup --samplesLabel ESC 2CLC --regionsLabel Boundary --plotHeight 6 --plotWidth 7 --colors blue red



computeMatrix reference-point -S Untreat.DI.NAto0.bw auxin-2days.DI.NAto0.bw washoff-2days.DI.NAto0.bw -R mESC_Dixon2012-raw_TADs.boundary.bed -b 600000 -a 600000 -o Boundary.CTCF.gz --referencePoint center --binSize 20000 -p 20

plotProfile -m Boundary.CTCF.gz -out Boundary.CTCF.pdf --perGroup --samplesLabel Untreat auxin2days washoff2days --regionsLabel Boundary --plotHeight 6 --plotWidth 7
