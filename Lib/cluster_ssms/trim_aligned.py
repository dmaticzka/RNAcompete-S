#!/home/hugheslab1/kate/bin/python

import sys
from numpy import *

alignedfile = sys.argv[1]
outfile = sys.argv[2]

cutoff = 0.5

totalseqs = 0

with open(alignedfile, 'r') as f:
	first_line = f.readline().rstrip()
	width = len(first_line)

#print 'width',width

numgap = [0]*width

af = open(alignedfile)
for line in af:
	kmer = array([bit for bit in line.rstrip()])
	is_gap = kmer[:] == '-'
	#print kmer,is_gap
	numgap = numgap + is_gap*1
	totalseqs = totalseqs + 1
af.close()

percgap = (1.0*numgap)/totalseqs

#print percgap

keepidx = percgap < cutoff

#print keepidx

startidx = 0
endidx = width

if keepidx[0] == False: # trim from the left
	for i in range(0,width):
		if keepidx[i] == False:
			startidx = i
		else:
			startidx = i
			break

if keepidx[width-1] == False: # trim from the right
	endidx = width-1
	for i in range(width-1,0,-1):
		#print i
		if keepidx[i] == False:
			endidx = i
		else:
			break

#print 'start',startidx,'endidx',endidx

af = open(alignedfile)
of = open(outfile,'w')
for line in af:
	kmer = line.rstrip()
	of.write(kmer[startidx:endidx]+'\n')
af.close()
of.close()
