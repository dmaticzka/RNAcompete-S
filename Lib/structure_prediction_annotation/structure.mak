RBP_FILE = $(shell head -2 ../info.tab | tail -1 | cut -f 1)
POOL_FILE = $(shell head -2 ../info.tab | tail -1 | cut -f 3)
RBPFILE = $(addprefix ../, $(RBP_FILE))
POOLFILE = $(addprefix ../, $(POOL_FILE))

RBP_centroid_structs.txt.gz: $(RBPFILE)
	zcat $< \
	| parallel --pipe -N1000 ../Lib/perl_utilities/fasta2tab.pl \
	| cut -f 1,2 \
	> ids.tmp;
	zcat $< \
	| ../Lib/perl_utilities/fasta2tab.pl \
	| cut -f 2 \
	| ../Lib/perl_utilities/paste.pl 'seq' - \
	| ../Lib/perl_utilities/tab2fasta.pl \
	| parallel --pipe -N10000 RNAfold -p --noPS \
	| ../Lib/structure_prediction_annotation/process_centroid_structure.pl \
	| sed -e 1d \
	| ../Lib/perl_utilities/fasta2tab.pl \
	| paste ids.tmp - \
	| cut -f 1,2,4 \
	| gzip > $@;
	rm -f seq_dp.ps seq_ss.ps ids.tmp;

RBP_secondary_structure.txt.gz: RBP_centroid_structs.txt.gz
	zcat $< \
	| cut -f 3 \
	> structs.tmp;
	\
	../Lib/structure_prediction_annotation/parse_secondary_structure structs.tmp annotated.tmp;
	\
	zcat $< \
	| cut -f 1,2 \
	| paste - annotated.tmp \
	| gzip > $@;
	rm -f structs.tmp annotated.tmp;

POOL_centroid_structs.txt.gz: $(POOLFILE)
	zcat $< \
        | ../Lib/perl_utilities/fasta2tab.pl \
        | cut -f 1,2 \
        > ids.tmp;
	zcat $< \
        | ../Lib/perl_utilities/fasta2tab.pl \
        | cut -f 2 \
        | ../Lib/perl_utilities/paste.pl 'seq' - \
        | ../Lib/perl_utilities/tab2fasta.pl \
        | parallel --pipe -N10000 RNAfold -p --noPS \
        | ../Lib/structure_prediction_annotation/process_centroid_structure.pl \
        | sed -e 1d \
        | ../Lib/perl_utilities/fasta2tab.pl \
        | paste ids.tmp - \
        | cut -f 1,2,4 \
        | gzip > $@;
	rm -f seq_dp.ps seq_ss.ps ids.tmp;

POOL_secondary_structure.txt.gz: POOL_centroid_structs.txt.gz
	zcat $< \
        | cut -f 3 \
        > structs.tmp;
	\
        ../Lib/structure_prediction_annotation/parse_secondary_structure structs.tmp annotated.tmp;
	\
        zcat $< \
        | cut -f 1,2 \
        | paste - annotated.tmp \
        | gzip > $@;
	rm -f structs.tmp annotated.tmp;
