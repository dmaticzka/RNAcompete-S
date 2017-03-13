#!/home/hugheslab1/kate/bin/python

from numpy import *
from scipy.sparse.linalg import eigs, ArpackNoConvergence
import sys



def compute_pfm_distance_and_alignment(kmer_seq_1,kmer_seq_2,kmer_struct_1,kmer_struct_2):
	a = array([bit for bit in kmer_seq_1])
	b = array([bit for bit in kmer_seq_2])
	a_struct = array([bit for bit in kmer_struct_1])
	b_struct = array([bit for bit in kmer_struct_2])
	
	#make sure a is longer
	swapped = False
	if len(a)<len(b):
		c=a.copy()
		a=b
		b=c
		
		c=a_struct.copy()
		a_struct=b_struct
		b_struct=c
		
		swapped=True #make sure a is longer
		
	#if swapped: print 'swapped!'
	
	np1=len(a)# no. of positions
	np2=len(b)
	dist=zeros(np1+np2)
	
	#print kmer1,kmer2
	#print a,b
	#print 'np1',np1,'np2',np2

	for i in arange(np1+np2):
		#print '==>start',i,j
		d=0
		if i<np2: # there is overlap on the left side
			#print 'left overlap',a[0:i],b[np2-i:np2]
			#print i,'num overlapping bases:',i,'num nonoverlapping bases',np2 - i + np1 - i
			#print i,'i<np2 overlap',sum(a[0:i]!=b[np2-i:np2]),d
			#print i,'i<np2 extra',0.75 * (np2 - i + np1 - i),d
			d+= sum(a[0:i]!=b[np2-i:np2]) + sum(a_struct[0:i]!=b_struct[np2-i:np2]) #add distance of overlapping bases from both sequence and structure kmers
			d+= 0.75 * (np2 - i + np1 - i) #add nonoverlapping distances
		elif i>np1: # there is overlap on the right side
			#print 'right overlap',a[i-np2:],b[0:np2+np1-i]
			#print i,'num overlapping bases:',np2+np1-i,'num nonoverlapping bases', i * 2 - np2 - np1
			#print i,'i>np1 overlap',sum(a[i-np2:] != b[0:np2+np1-i]),d
			#print i,'i>np1 extra',0.75 * (i * 2 - np2 - np1),d
			d+= sum(a[i-np2:] != b[0:np2+np1-i]) + sum(a_struct[i-np2:] != b_struct[0:np2+np1-i])
			d+= 0.75 * (i * 2 - np2 - np1)
		else: # no overlap, smaller kmer is completely within larger kmer
			#print 'completely contained'
			#print i,np1,np2,a[i-np2:i],b[:]
			#print i,'i>=np2 overlap',sum(a[i-np2:i]!=b[:]),d
			#print i,'i>=np2 extra',0.75 * (np1 - np2),d
			d+=sum(a[i-np2:i]!=b[:]) + sum(a_struct[i-np2:i]!=b_struct[:])
			d+= 0.75 * (np1 - np2) # extra bases = length that longer kmer exceeds shorter kmer
		dist[i]=d
		#print i,'final d',d
	align = argmin(dist)
# 	print(dist)
# 	print align
# 	print np1,np2
# 	if(swapped):
# 		print '-'*(align-np2) + kmer_seq_1 + '-'*(np1 - align)
# 		print '-'*(np2-align) + kmer_seq_2 + '-'*(align-np1)
# 		print '-'*(align-np2) + kmer_struct_1 + '-'*(np1 - align)
# 		print '-'*(np2-align) + kmer_struct_2 + '-'*(align-np1)
# 	else:
# 		print '-'*(np2-align) + kmer_seq_1 + '-'*(align-np1)
# 		print '-'*(align-np2) + kmer_seq_2 + '-'*(np1 - align)
# 		print '-'*(np2-align) + kmer_struct_1 + '-'*(align-np1)
# 		print '-'*(align-np2) + kmer_struct_2 + '-'*(np1 - align)
	return min(dist),align


def align_kmers(kmers_seq,kmers_struct):
	#print(kmers)
	m = kmers_seq.shape[0] # no. of kmers

	# calculate distance matrix
	D = zeros((m,m))
	A = zeros((m,m))
	for i in arange(m):
		for j in arange(m): #arange(i+1,m):
			dist,alignment = compute_pfm_distance_and_alignment(kmers_seq[i],kmers_seq[j],kmers_struct[i],kmers_struct[j])
			D[i,j]= dist
			#D[j,i]=D[i,j]
			A[i,j]= alignment
			#A[j,i]=A[i,j]
	
	#print kmers[0],A[0,:]
	
	# find kmer with the lowest total distance (summed across all other kmers)
	seedkmer = 0
	seedkmertotaldist = inf
	for i in arange(m):
		totaldist = 0
		for j in arange(m):
			if i == j:
				continue
			totaldist += D[i,j]
		#print i,kmers[i],totaldist,seedkmer,seedkmertotaldist
		if totaldist < seedkmertotaldist:
			seedkmer = i
			seedkmertotaldist = totaldist
	
	# find the left padding on the seed kmer, we will add this to the padding
	# for all other kmers
	seedpadding = 0
	for i in arange(m):
		curpadding = 0
		np1 = len(kmers_seq[seedkmer])
		np2 = len(kmers_seq[i])
		#print 'kmer1',kmers_seq[seedkmer],'kmer2',kmers_seq[i],'np1',np1,'np2',np2
		if np1 > np2:
			curpadding = np1-A[seedkmer,i]
		else:
			curpadding = abs(A[seedkmer,i]-np1)
		#print 'curpadding:',curpadding
		if curpadding > seedpadding:
			seedpadding = curpadding
	
	
	#print 'seedkmer',kmers_seq[seedkmer],kmers_struct[seedkmer]
	#print 'seedpadding',seedpadding
	
	# get left padding for all k-mers using alignment to seed kmer
	leftpadding = [0]*m
	np1 = len(kmers_seq[seedkmer])
	for i in arange(m):
		np2 = len(kmers_seq[i])
		align = A[seedkmer,i]
		if i == seedkmer:
			leftpadding[i] = seedpadding
		elif len(kmers_seq[seedkmer]) < len(kmers_seq[i]): #swapped
			leftpadding[i] = np1-align + seedpadding
		else:
			leftpadding[i] = align-np1 + seedpadding
	
	# figure out right padding
	maxlen = 0
	for i in arange(m):
		curlen = leftpadding[i] + len(kmers_seq[i])
		if curlen > maxlen:
			maxlen = curlen
	
	rightpadding = [0]*m	
	for i in arange(m):
		rightpadding[i] = maxlen - (leftpadding[i] + len(kmers_seq[i]))
	
	aligned_seq = []
	aligned_struct = []
	
	for i in arange(m):
		aligned_seq.append('-'*(leftpadding[i]) + kmers_seq[i] + '-'*(rightpadding[i]))
		aligned_struct.append('-'*(leftpadding[i]) + kmers_struct[i] + '-'*(rightpadding[i]))
		#print '-'*(leftpadding[i]) + kmers_seq[i] + '-'*(rightpadding[i]) + '\t' + '-'*(leftpadding[i]) + kmers_struct[i] + '-'*(rightpadding[i])
	
	return aligned_seq,aligned_struct



if __name__ == "__main__":
	#print compute_pfm_distance_and_alignment('AGA','AAAGAG','HHH','RRHHHH')
	kmers_seq=loadtxt(sys.argv[1],dtype='string')
	kmers_struct=loadtxt(sys.argv[2],dtype='string')
	aligned_seq,aligned_struct= align_kmers(kmers_seq,kmers_struct)
	
	of = open(sys.argv[3], 'w')
	of.write('\n'.join(aligned_seq))
	of.close()
	
	of = open(sys.argv[4], 'w')
	of.write('\n'.join(aligned_struct))
	of.close()

	

