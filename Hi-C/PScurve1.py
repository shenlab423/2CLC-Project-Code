import sys
import numpy as np

chroms = {}
bininfo = {}
with open(sys.argv[1]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        chroms[line_s[-1]] = line_s[0]
        info1 = '-'.join(line_s[0:3])
        bininfo[info1] = int(line_s[-1])

alls = 0
distance_x = {}
interaction = {}
with open(sys.argv[2]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = '_'.join(line_s[0:2])
        locusinfo2 = '_'.join(line_s[1::-1])
        interaction[locusinfo] = line_s[2]
        interaction[locusinfo2] = line_s[2]

        froms = line_s[0]
        tos = line_s[1]
        if chroms[froms] == chroms[tos] and chroms[froms]=="chr1":
            alls += float(line_s[2])
            chrom = chroms[froms]
            distance = abs(int(froms) - int(tos))

            if distance in distance_x:
                distance_x[distance].append(float(line_s[2]))
            else:
                distance_x[distance] = [float(line_s[2])]

distance_x2 = {}
distance_x3 = {}
for distance,value2 in distance_x.items():
    distance_x2[distance] = np.mean(value2)
    distance_x3[distance] = np.sum(value2)

out2 = [(k,distance_x2[k]) for k in sorted(distance_x2.keys())]

for x,y in out2:
    print("%s\t%s\t%s\t%s"%(x,y,distance_x3[x],alls)) 


