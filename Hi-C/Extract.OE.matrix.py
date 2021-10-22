import sys
import numpy as np
#chrom start end
#start and end must be binsize*int

chroms = {}
bininfo = {}
with open(sys.argv[1]) as infile2:
	for line in infile2:
		line = line.strip()
		line_s = line.split('\t')
		chroms[line_s[-1]] = line_s[0]
		info1 = '-'.join(line_s[0:3])
		bininfo[info1] = int(line_s[-1])

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
        if chroms[froms] == chroms[tos]:
            chrom = chroms[froms]
            distance = abs(int(froms) - int(tos))
            if not chrom in distance_x:
                distance_x[chrom] = {}

            if distance in distance_x[chrom]:
                distance_x[chrom][distance].append(float(line_s[2]))
            else:
                distance_x[chrom][distance] = [float(line_s[2])]

distance_x2 = {}
for chrom,value in distance_x.items():
    distance_x2[chrom] = {}
    for distance,value2 in value.items():
        distance_x2[chrom][distance] = np.mean(value2)

distance_y = {}
interaction2 = {}
with open(sys.argv[3]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = '_'.join(line_s[0:2])
        locusinfo2 = '_'.join(line_s[1::-1])
        interaction2[locusinfo] = line_s[2]
        interaction2[locusinfo2] = line_s[2]

        froms = line_s[0]
        tos = line_s[1]
        if chroms[froms] == chroms[tos]:
            chrom = chroms[froms]
            distance = abs(int(froms) - int(tos))
            if not chrom in distance_y:
                distance_y[chrom] = {}

            if distance in distance_y[chrom]:
                distance_y[chrom][distance].append(float(line_s[2]))
            else:
                distance_y[chrom][distance] = [float(line_s[2])]

distance2 = {}
for chrom,value in distance_y.items():
    distance2[chrom] = {}
    for distance,value2 in value.items():
        distance2[chrom][distance] = np.mean(value2)

		
#print(interaction2["53703_53727"])
#print(interaction["110185_53570"])
#matrix = np.zeros((1,50,50))

binsize = 20000

loucs = sys.argv[4]

chrom_start_info = loucs.split(':')[0] + '-' + loucs.split(':')[1].split('-')[0] + '-' + str(int(loucs.split(':')[1].split('-')[0])+binsize)
chrom_end_info = loucs.split(':')[0] + '-' +  str(int(loucs.split(':')[1].split('-')[1])-binsize) + '-' + loucs.split(':')[1].split('-')[1]
bins_before = int(bininfo[chrom_start_info])
bins_after = int(bininfo[chrom_end_info]) + 1
lengths = bins_after - bins_before
matrix = np.full((lengths , lengths ), np.nan)
matrix2 = np.full((lengths , lengths ), np.nan)
for i in range(bins_before,bins_after):
	for j in range(bins_before,bins_after):
		distance = abs(i - j)
		index_i = i - bins_before
		index_j = j - bins_before
		strs = str(i) + '_' + str(j)
		#print(index_i,index_j)
		#print(float(interaction[strs]))
		if index_j > index_i:
			matrix[index_i,index_j] = np.nan
		else:
			#if abs(index_j - index_i) > (1 / 6 * lengths):
			#	matrix[index_i,index_j] = np.nan
			#else:
			if strs in interaction:
				matrix[index_i,index_j] = float(interaction[strs]) / float(distance_x2[loucs.split(':')[0]][distance])
for i in range(bins_before,bins_after):
	for j in range(bins_before,bins_after):
		distance = abs(i - j)
		index_i = i - bins_before
		index_j = j - bins_before
		strs = str(i) + '_' + str(j)
		#print(index_i,index_j)
		#print(float(interaction[strs]))
		if index_j > index_i:
			matrix2[index_i,index_j] = np.nan
		else:
			#if abs(index_j - index_i) > (1 / 6 * lengths):
			#	matrix2[index_i,index_j] = np.nan
			#else:
			if strs in interaction2:
				matrix2[index_i,index_j] = float(interaction2[strs])  / float(distance2[loucs.split(':')[0]][distance])
				
np.savetxt(sys.argv[5]+'.txt',matrix)
np.savetxt(sys.argv[6]+'.txt',matrix2)


