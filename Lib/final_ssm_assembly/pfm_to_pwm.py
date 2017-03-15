import os,sys,csv

from pfm_util import *

infile = sys.argv[1]
outfile = sys.argv[2]
num_sites = float(sys.argv[3])


pfm = read_pfm(infile)

pwm = pfm_to_pwm(pfm,num_sites)

write_pfm(pwm,outfile)

