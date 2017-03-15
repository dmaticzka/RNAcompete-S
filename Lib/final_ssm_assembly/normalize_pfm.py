import os,sys,csv

from pfm_util import *

infile = sys.argv[1]
outfile = sys.argv[2]


pfm = read_pfm(infile)

pfm_normalized = norm_pfm(pfm)

write_pfm(pfm_normalized,outfile)

