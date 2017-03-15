
import os,sys,csv

import pfm_util


def pfm_combine(pfm_1,pfm_2,offset):
	# set up output PFM
	pfm_final = {}
	alphabet = pfm_1.keys()
	pfm_length_1 = len(pfm_1[alphabet[0]])
	pfm_length_2 = len(pfm_2[alphabet[0]])

	if offset <= 0 and pfm_length_1 <= pfm_length_2:
		pfm_length_final = abs(offset)+max(pfm_length_1,pfm_length_2)
	elif offset <=0 and pfm_length_1 >= pfm_length_2:
		pfm_length_final = abs(offset)+max(pfm_length_1,pfm_length_2)
	elif offset >=0 and pfm_length_1 <= pfm_length_2:
		pfm_length_final = abs(offset)+max(pfm_length_1,pfm_length_2)
	else: #offset >=0 and pfm_length_1 >= pfm_length_2
		pfm_length_final = abs(offset)+max(pfm_length_1,pfm_length_2)


	for base in alphabet:
		pfm_final[base] = [None] * (pfm_length_final)


	# print('pfm_lengths',pfm_length_1,pfm_length_2)
	# print('pfm_length_final',pfm_length_final)
	# print(pfm_final)

	max_pos_pfm_1 = 0 # highest position that we've dealt with
	max_pos_pfm_2 = 0
	max_pos_final = 0

	# add positions that are left of the overlapping region

	if offset > 0: # extra is from PFM 1
		for pos in range(0, offset):
			# print('adding left pos',pos)
			for base in alphabet:
				pfm_final[base][pos] = pfm_1[base][pos]
			max_pos_pfm_1 = pos
			max_pos_final = pos
	elif offset < 0: # extra is from PFM 2
		for pos in range(0, -offset):
			# print('adding left pos',pos)
			for base in alphabet:
				pfm_final[base][pos] = pfm_2[base][pos]
			max_pos_pfm_2 = pos
			max_pos_final = pos

	# print('after left pos')
	# print(pfm_final)

	# figure out borders of overlapping region
	# add overlapping positions

	# print('adding overlap positions')

	offset_1 = 0
	offset_2 = 0
	offset_final = 0

	if offset > 0:
		# positions in PFM 1 are +offset
		# positions in PFM 2 are right
		# positions in final are +offset
		# start position in PFM 1 is offset
		# start position in PFM 2 is 0
		for pos in range(0,min(pfm_length_1-offset,pfm_length_2)):
			for base in alphabet:
				print('a',base,pos,pos+offset)
				print pfm_1
				print pfm_2
				pfm_final[base][pos+offset] = pfm_1[base][pos+offset] + pfm_2[base][pos]
			max_pos_pfm_2 = pos
			max_pos_pfm_1 = pos+offset
			max_pos_final = pos+offset
			offset_1 = offset
			offset_2 = 0
			offset_final = offset
	elif offset < 0:
		# positions in PFM 1 are + -offset
		# positions in PFM 2 are right
		# start position in PFM 1 is 0
		# start position in PFM 2 is -offset
		for pos in range(0,min(pfm_length_1,pfm_length_2+offset)):
			for base in alphabet:
				# print('b',base,pos,pos-offset)
				pfm_final[base][pos-offset] = pfm_1[base][pos] + pfm_2[base][pos-offset]
			max_pos_pfm_2 = pos-offset
			max_pos_pfm_1 = pos
			max_pos_final = pos-offset
			offset_1 = 0
			offset_2 = -offset
			offset_final = -offset
	else:
		# offset is zero
		# positions in PFM 1 and 2 are right
		# start position in both are 0
		for pos in range(0,min(pfm_length_1,pfm_length_2)):
			for base in alphabet:
				# print('x',base,pos)
				pfm_final[base][pos] = pfm_1[base][pos] + pfm_2[base][pos]
			max_pos_pfm_2 = pos
			max_pos_pfm_1 = pos
			max_pos_final = pos
			offset_1 = 0
			offset_2 = 0
			offset_final = 0

	# print('after overlap')
	# print(pfm_final)
	# print('max_pos_1',max_pos_pfm_1)
	# print('length_1',pfm_length_1)
	# print('max_pos_2',max_pos_pfm_2)
	# print('length_2',pfm_length_2)
	# print('max_pos_final',max_pos_final)
	# print('')
	# print('offset_1',offset_1)
	# print('offset_2',offset_2)
	# print('offset_final',offset_final)
	# print('')

	# add positions that are right of the overlapping region

	# print('adding right pos')

	if max_pos_pfm_1 < pfm_length_1-1:
		for pos in range(max_pos_pfm_1+1,pfm_length_1):
			for base in alphabet:
				# print('c',base,pos)
				pfm_final[base][pos+offset_final-offset_1] = pfm_1[base][pos]
				max_pos_pfm_1 = pos
				max_pos_final = pos+abs(offset)
	elif max_pos_pfm_2 < pfm_length_2-1:
		for pos in range(max_pos_pfm_2+1,pfm_length_2):	
			for base in alphabet:
				# print('d',base,pos)
				pfm_final[base][pos+offset_final-offset_2] = pfm_2[base][pos]
				max_pos_pfm_2 = pos
				max_pos_final = pos+abs(offset)

	return pfm_final

# infile1 = sys.argv[1]
# infile2 = sys.argv[2]
# outfile = sys.argv[3]
# offset = int(sys.argv[4])

# Offset is difference between start positions and is relative to pfm 1, so like this:

# [X X X X] pfm 1
#   [X X X X X] pfm 2
# offset = 1

# [X X X X] pfm 1
# [X X X] pfm 2
# offset = 0

#     [X X X X] pfm 1
# [X X X X X] pfm 2
# offset = -2

# read in PFMs
# pfm_1 = read_pfm(infile1)
# pfm_2 = read_pfm(infile2)
# 
# pfm_final = pfm_combine(pfm_1,pfm_2,offset)
# 
# write_pfm(pfm_final,outfile)