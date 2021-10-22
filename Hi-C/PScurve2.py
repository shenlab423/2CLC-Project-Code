import sys
import numpy as np
from pylab import *
from scipy import optimize
from numpy import log10
from sklearn import linear_model
from math import *
from math import pow
from random import shuffle
from collections import Counter


dict1 = {}
with open(sys.argv[1]) as infile:
    for line in infile:
        line = line.strip()
        splits = line.split("\t")
        distance = int(splits[0]) * 5000
        dict1[distance] = float(splits[2])

dict2 = {}
#log_bins = [0,20000,100000,600000,3600000,20000000,115000000]
#log_bins = [0,20000,40000,80000,160000,320000,640000,1280000,2560000,5120000,10240000,20480000,40960000,81920000,163840000]
log_bins = [4000,8000,17000,34000,68000,136000,271000,543000,1050000,2100000,4190000,8390000,16780000,33550000,67110000,134220000]

for i in range(len(log_bins)):
    y = i + 1
    if y < len(log_bins):
        start = log_bins[i]
        end = log_bins[y]
        sums = 0.00000
        nums = 0
        for distance,value in dict1.items():
            if distance >= start and distance < end:
                sums += value
                nums += 1
        dict2[start] = sums
    else:
        start = log_bins[i]
        end = 0
        sums = 0.00000
        nums = 0
        for distance,value in dict1.items():
            if distance >= start:
                sums += value
                nums += 1
        dict2[start] = sums

    print("%s\t%s\t%s"%(start,end,sums))
