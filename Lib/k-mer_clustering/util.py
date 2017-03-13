from numpy import *
from scipy.sparse.linalg import eigs, ArpackNoConvergence

# computes distance between two kmers by counting the number 
# of matching bases and then converting it into a distance metric
#
# kmer1: input kmer1
# kmer2: input kmer2
# return distance score
def compute_pfm_distance(kmer1,kmer2):
	a = array([bit for bit in kmer1])
	b = array([bit for bit in kmer2])
	if len(a)<len(b):c=a.copy();a=b;b=c;		
	np1=len(a)# no. of positions
	np2=len(b)
	kern=zeros(np1+np2)

	j=0		
	for i in arange(np1+np2-1):
		d=0
		if i<np2:
			d+= sum(a[0:i+1]==b[np2-i-1:np2])
		else:
			j+=1
			d+=sum(a[j:j+min(np2,np1-j)]==b[0:min(np2,np1-j)])

		kern[i]=d	
	return sqrt(np1+np2-2*max(kern))

	

# computes distance between two kmers (by comparing the 0/1 PFMs)
# kmer1: input kmer1
# kmer2: input kmer2
# return distance score
def compute_pfm_distance1(kmer1,kmer2):
	a = array([bit for bit in kmer1])
	b = array([bit for bit in kmer2])
	if len(a)<len(b):c=a.copy();a=b;b=c;		
	np1=len(a)# no. of positions
	np2=len(b)
	dist=zeros(np1+np2)

	j=0		
	for i in arange(np1+np2-1):
		d=0
		if i<np2:
			d+= sum(a[0:i+1]!=b[np2-i-1:np2])
			d+= 0.75*(np1-(i+1))#len(a[i+1::])
			if np2-i-1>0: d+=0.75*(np2-i-1)#len(b[0:np2-1-1])
		else:
			j+=1
			d+=sum(a[j:j+min(np2,np1-j)]!=b[0:min(np2,np1-j)])
			d+=0.75*j#len(a[0:j])
			if j>np1-np2: 
				d+=0.75*(np2-(np1-j))#0.75*len(b[np1-j::])
			else:
				d+=0.75*(np1-(np2+j))#0.75*len(a[np2+j::])
		dist[i]=d
	dist[np1+np2-1] = 0.75*(np1+np2)
	
	return min(dist)



# check if the distance matrix D results in a PSD similarity matrix K
# if not transform the distance matrix, and compute its corresponding K
def compute_kernel(D):
	m = D.shape[0]
	Q = eye(m)-ones((m,m))*1./m
	K = -0.5*dot(dot(Q,D),Q)
	fail = 0

	try:
		eval,evec=eigs(K,k=1,which='SM') # smallest eigenvalue
		eval=eval.real
	except ArpackNoConvergence:
		eval=0;fail=1

	if eval < 0: # K is not PSD
		D = D - 2*eval*(ones((m,m))-eye(m))		# transform D following the PAMI paper
		K = -0.5*dot(dot(Q,D),Q)	
	
	return D, K, fail