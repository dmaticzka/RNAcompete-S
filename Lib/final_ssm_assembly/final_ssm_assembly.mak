DISTFILE = all_dist_modes.txt
PFMDIR = .

targets = 

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp) $(wildcard *.cache)

a:
	echo $(DISTFILE)

input.txt: $(DISTFILE)
	cat $< \
	| sed 's/\t/_/' \
	| ../Lib/perl_utilities/row_stats.pl -abs -mean \
	| paste - $(DISTFILE) \
	| sed -e 1d \
	| ../Lib/perl_utilities/cut.pl -f 3,4,5,2 \
	| sort -nk 4 -k 1 \
	| tr : '\t' \
	| ../Lib/perl_utilities/cut.pl -f 2,4,5,6 \
	| ../Lib/perl_utilities/cap.pl key,other,distmode,absdistmode \
	> $@;

pfm_assembled_seq.txt: input.txt
	python ../Lib/final_ssm_assembly/assemble_pfm.py $< $@ $(PFMDIR) pfm_seq.txt ;

pfm_assembled_struct.txt: input.txt
	python ../Lib/final_ssm_assembly/assemble_pfm.py $< $@ $(PFMDIR) pfm_struct.txt ;

pfm_assembled_seq_normalized.txt: pfm_assembled_seq.txt
	python ../Lib/final_ssm_assembly/normalize_pfm.py $< $@; \
	rm -f pfm_assembled_seq.txt;

pfm_assembled_struct_normalized.txt: pfm_assembled_struct.txt
	python ../Lib/final_ssm_assembly/normalize_pfm.py $< $@; \
	rm -f pfm_assembled_struct.txt;

pwm_seq.txt: pfm_assembled_seq_normalized.txt
	python ../Lib/final_ssm_assembly/pfm_to_pwm.py $< $@ 100;

pwm_struct.txt: pfm_assembled_struct_normalized.txt
	python ../Lib/final_ssm_assembly/pfm_to_pwm.py $< $@ 100;

