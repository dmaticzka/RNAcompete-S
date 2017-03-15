import sys,csv
from pfm_util import *
import combine_pfms

infile = sys.argv[1]
outfile = sys.argv[2]
pfmdir = sys.argv[3]
pfmfilename = sys.argv[4]
# $(HOME)/random_pulldown_analysis_new/cluster_pfms_and_logos_withweights/$(RBP_ID)/

clust_ids_done = []
offsets = {}
cur_pfm = None

with open(infile) as f:
	reader = csv.reader(f, delimiter='\t')
	headerline = reader.next()
	for row in reader:
		key = row[0]
		other = row[1]
		distmode = int(row[2])
		if key in clust_ids_done and other in clust_ids_done:
			#print('skipping')
			continue
		if cur_pfm is None:		
			print('key',key,'other',other,'distmode',distmode)
			pfm_1 = read_pfm(pfmdir + '/' + key + '/' + pfmfilename)
			pfm_2 = read_pfm(pfmdir + '/' + other + '/' + pfmfilename)
			cur_pfm = combine_pfms.pfm_combine(pfm_1,pfm_2,distmode)
			if distmode > 0:
				offsets[key] = 0
				offsets[other] = distmode
			else:
				offsets[key] = -distmode
				offsets[other] = 0
			clust_ids_done.append(key)
			clust_ids_done.append(other)
		else:
			if other in clust_ids_done:
				#print('skipping')
				continue			
			print('key',key,'other',other,'distmode',distmode)
			# key is in clust_ids_done
			addoffset = offsets[key]
			print('addoffset',addoffset,'totaldist',distmode+addoffset)
			pfm_1 = cur_pfm
			pfm_2 = read_pfm(pfmdir + '/' + other + '/' + pfmfilename)
			cur_pfm = combine_pfms.pfm_combine(pfm_1,pfm_2,distmode+addoffset)
			
			if distmode+addoffset > 0:
				offsets[key] = 0
				offsets[other] = distmode+addoffset
			else:
				for k in offsets.keys():
					offsets[k] = offsets[k] - (distmode+addoffset)
				offsets[other] = 0
			
			offsets[other] = distmode+addoffset
			clust_ids_done.append(other)			
		#write_pfm(cur_pfm,'testout'+key+other+'.txt')


f.close()

write_pfm(cur_pfm,outfile)