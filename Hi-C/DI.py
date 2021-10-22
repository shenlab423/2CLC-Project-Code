import sys
import math
import numpy as np

chrominfo = {}
chrombins = {}
chroms = {}

with open(sys.argv[1]) as infile:
    for line in infile:
        line = line.strip()
        line_s = line.split('\t')
        chroms[line_s[-1]] = line_s[0]
        locusinfo = line_s[0]+'\t'+line_s[1]+'\t'+line_s[2]
        if line_s[0] in chrominfo:
            chrominfo[line_s[0]][int(line_s[3])] = locusinfo
            chrombins[line_s[0]].append(int(line_s[3]))
        else:
            chrominfo[line_s[0]] = {}
            chrominfo[line_s[0]][int(line_s[3])] = locusinfo
            chrombins[line_s[0]] = []
            chrombins[line_s[0]].append(int(line_s[3]))


interaction = {}
with open(sys.argv[2]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = '_'.join(line_s[0:2])
        locusinfo2 = '_'.join(line_s[1::-1])
        interaction[locusinfo] = float(line_s[2])
        interaction[locusinfo2] = float(line_s[2])



for key,value_list in chrombins.items():
    with open(key+".DI.bedGraph",'w') as output:
        strs = '\t'

        maxvalue = max(value_list)
        for m in value_list:
            strs = chrominfo[key][m]
            if m <= 50 or (m + 50) > maxvalue:
                output.write(strs.strip()+'\tNA\n')
                continue
            start = m - 50
            end = m + 51
            
            A = 0
            for x in range(start,m):
                y = m
                distance_t = abs(x - y)
                if str(x)+'_'+str(y) in interaction:
                    A += interaction[str(x)+'_'+str(y)]

            B = 0
            for x in range(m+1,end):
                y = m
                distance_t = abs(x - y)
                if str(x)+'_'+str(y) in interaction:
                    B += interaction[str(x)+'_'+str(y)]

            if A + B > 0:
                E = 0.5 * (A + B )
                DI = ((B - A)/abs(B-A)) * ((A-E) * (A-E) + (B-E) * (B-E)) / E
            else:
                DI = "NA"
            output.write(strs.strip()+'\t'+str(DI)+'\n')
