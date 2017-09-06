#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift @ARGV;
my $label = shift @ARGV;
my $minlen = shift @ARGV || 2;
my $maxlen = shift @ARGV || 8;

my $seqNamespace = "";
my $structNamespace = "";
my $jointNamespace = "";


open(INSTREAM, $file) or die "couldn't open\n";


while (<INSTREAM>){
    chomp; chop while /\r/;
        my ($id, $seq, $struct) = split("\t");
        #my ($id, $seq, $struct) = split("\t");
        #my $label = 1;

        my @seqChars = split (//, $seq);
        my @structChars = split (//, $struct);

        for (my $len=$minlen; $len<$maxlen+1; $len++) {

            if ($len == 2){
                $jointNamespace = "q";
            } elsif ($len == 3) {
                $jointNamespace = "r";
            } elsif ($len == 4) {
                $jointNamespace = "q";
            } elsif ($len == 5) {
                $jointNamespace = "s";
            } elsif ($len == 6) {
                $jointNamespace = "u";
            } elsif ($len == 7) {
                $jointNamespace = "v";
            } elsif ($len == 8) {
                $jointNamespace = "w";
            } else {
                $jointNamespace = "x";
            }

            my @joint = ();

            for (my $i=0; $i<scalar @seqChars-$len+1; $i++) {
                my $curSeq = substr ($seq, $i, $len);
                my $curStruct = substr ($struct, $i, $len);
                push(@joint,"${curSeq}_$curStruct");
            }
            if ($len == $minlen) {
                print $label . " " . $id;
            }
            print " |$jointNamespace " . join(" ",@joint) . " ";
            if ($len == $maxlen) {
                print "\n"
            }
        }
}

exit(0);
