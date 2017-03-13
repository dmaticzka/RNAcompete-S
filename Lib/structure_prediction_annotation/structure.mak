RBP_FILE = $(shell head -2 ../info.tab | tail -1 | cut -f 1)
POOL_FILE = $(shell head -2 ../info.tab | tail -1 | cut -f 3)
RBPFILE = $(addprefix ../, $(RBP_FILE))
POOLFILE = $(addprefix ../, $(POOL_FILE))

RBP_centroid_structs.txt: $(RBPFILE)
	zcat $< \
	| ../Lib/perl/fasta2tab.pl \
	| cut -f 1,2 \
	> ids.tmp;
	zcat $< \
	| ../Lib/perl/fasta2tab.pl \
	| cut -f 2 \
	| ../Lib/perl/paste.pl 'seq' - \
	| ../Lib/perl/tab2fasta.pl \
	| RNAfold -p --noPS \
	| ../Lib/structure_prediction_annotation/process_centroid_structure.pl \
	| sed -e 1d \
	| ../Lib/perl/fasta2tab.pl \
	| paste ids.tmp - \
	| cut -f 1,2,4 \
	> $@;
	rm -f seq_dp.ps seq_ss.ps ids.tmp;

RBP_centroid_structs.txt.gz: RBP_centroid_structs.txt
	gzip $<;

RBP_secondary_structure.txt: RBP_centroid_structs.txt.gz
	zcat $< \
	| cut -f 3 \
	> structs.tmp;
	\
	../Lib/structure_prediction_annotation/parse_secondary_structure_v2 structs.tmp annotated.tmp;
	\
	zcat $< \
	| cut -f 1,2 \
	| paste - annotated.tmp \
	> $@;
	rm -f structs.tmp annotated.tmp;

RBP_secondary_structure.txt.gz: RBP_secondary_structure.txt
	gzip $<;

POOL_centroid_structs.txt: $(POOLFILE)
	zcat $< \
        | ../Lib/perl/fasta2tab.pl \
        | cut -f 1,2 \
        > ids.tmp;
	zcat $< \
        | ../Lib/perl/fasta2tab.pl \
        | cut -f 2 \
        | ../Lib/perl/paste.pl 'seq' - \
        | ../Lib/perl/tab2fasta.pl \
        | RNAfold -p --noPS \
        | ../Lib/structure_prediction_annotation/process_centroid_structure.pl \
        | sed -e 1d \
        | ../Lib/perl/fasta2tab.pl \
        | paste ids.tmp - \
        | cut -f 1,2,4 \
        > $@;
	rm -f seq_dp.ps seq_ss.ps ids.tmp;

POOL_centroid_structs.txt.gz: POOL_centroid_structs.txt
	gzip $<;

POOL_secondary_structure.txt: POOL_centroid_structs.txt.gz
	zcat $< \
        | cut -f 3 \
        > structs.tmp;
	\
        ../Lib/structure_prediction_annotation/parse_secondary_structure_v2 structs.tmp annotated.tmp;
	\
        zcat $< \
        | cut -f 1,2 \
        | paste - annotated.tmp \
        > $@;
	rm -f structs.tmp annotated.tmp;

POOL_secondary_structure.txt.gz: POOL_secondary_structure.txt
	gzip $<;
