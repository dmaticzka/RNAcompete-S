#!/home/hugheslab1/kate/bin/python


# reads two tab-delimited files and sums the contents. Assumes they are the same size
# (number of rows and columns)

import sys
import csv
from itertools import izip
import numpy as np

infile1 = sys.argv[1]
infile2 = sys.argv[2]

outfile = sys.argv[3]


with open(infile1) as if1: 
	with open(infile2) as if2: 
		with open(outfile, 'w') as of:
			reader1 = csv.reader(if1,delimiter='\t')
			reader2 = csv.reader(if2,delimiter='\t')
			for x, y in izip(reader1, reader2):
				z = np.add(np.array(map(float,x)),np.array(map(float,y)))
				of.write("\t".join(map(str,z)))
				of.write("\n")
