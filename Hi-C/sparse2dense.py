import sys

chrominfo = {}
chrombins = {}

with open(sys.argv[1]) as infile:
    for line in infile:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = line_s[3]+'|mm9|'+line_s[0]+':'+line_s[1]+'-'+line_s[2]
        if line_s[0] in chrominfo:
            chrominfo[line_s[0]][line_s[3]] = locusinfo
            chrombins[line_s[0]].append(line_s[3])
        else:
            chrominfo[line_s[0]] = {}
            chrominfo[line_s[0]][line_s[3]] = locusinfo
            chrombins[line_s[0]] = []
            chrombins[line_s[0]].append(line_s[3])

interaction = {}
with open(sys.argv[2]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = '_'.join(line_s[0:2])
        locusinfo2 = '_'.join(line_s[1::-1])
        interaction[locusinfo] = line_s[2]
        interaction[locusinfo2] = line_s[2]

for key,value_list in chrombins.items():
    with open(key+".dense.matrix.txt",'w') as output:
        strs = '\t'
        for m in value_list:
            strs = strs + chrominfo[key][m] + '\t'
        strs = '\t' + strs.strip()
        output.write(strs+'\n')
        
        for m in value_list:
            strs = chrominfo[key][m] + '\t'
            for n in value_list:
                if m+'_'+n in interaction:
                    strs = strs + interaction[m+'_'+n] + '\t'
                else:
                    strs = strs + 'nan\t'
            strs = strs.strip()
            output.write(strs+'\n')
