#!/usr/bin/perl


use strict;
require "libfile.pl";
use List::Util qw(sum);

$| = 1;

my @flags   = (
                  [    '-s', 'scalar',     0, 1]
                 ,[    '-c', 'scalar',     0, 1]
                 ,[    '-g', 'scalar',     0, 1]
                 ,['--file', 'scalar',   '-', undef]
              );

my %args = %{&parseArgs(\@ARGV, \@flags)};

if(exists($args{'--help'}))
{
   print STDOUT <DATA>;
   exit(0);
}

my $file  = $args{'--file'};
my $bSmallSampleCor = $args{'-s'};
my $bGapAdjust = !$args{'-g'};
my $bGetCounts = $args{'-c'};

die "error: can't use small sample correction with counts!" if ($bSmallSampleCor && $bGetCounts);

#open(my $fh, $file) or die "could not open file '$file'";
my $rarKmersZscores = readFileColumns($file,[0,1]);

my $nSeqs = scalar @{$rarKmersZscores};
my @aKmers;
my @aZscores;

for(my $i=0; $i<$nSeqs; $i++){
	push(@aKmers,$rarKmersZscores->[$i]->[0]);
	push(@aZscores,$rarKmersZscores->[$i]->[1]);
}

#warn @aKmers,"\n";
#warn @aZscores,"\n";

my $avgZ = sum(@aZscores)/(scalar @aZscores);
#warn $avgZ,"\n";

#my $sAlphabetStr = `cat  $file | cut -f 1 | tr -d '\n' | transpose.pl -d ''  -q | sort -u | grep -v '-' | transpose.pl -q`;
#chomp $sAlphabetStr;
#my @aAlphabet = split("\t",$sAlphabetStr);

my @aAlphabet = ('A','C','G','U');

#warn join(';',@aAlphabet),"\n";
#warn join(';',@aZscores),"\n";


my $rPFM = get_pfm_from_aligned_kmers(\@aKmers,\@aZscores,$bSmallSampleCor,$avgZ);

my $output = '';
foreach my $position (@{$rPFM}){
	foreach my $b (@aAlphabet){
		$output .= $position->{$b} . "\t";
	}
	chop $output;
	$output .= "\n";
}
print $output;

# ===============
# SUBROUTINES

sub get_pfm_from_aligned_kmers{
	my ($raKmers,$raZscores,$bSmallSampleCor,$avgZ) = @_;
	#make frequency matrix
	my $nWidth = length($raKmers->[0]); #assume alignments are all same length
	my $rPFM;
	my @anTotals = (0) x $nWidth;
	# count occurrence of each base at each position
	for(my $i=0; $i<$nWidth; $i++){ #position in width
		foreach my $b (@aAlphabet){
			#warn $i,$b,"\n";
			$rPFM->[$i]->{$b} = 0;
			my $n = 0;
			foreach my $sKmer (@{$raKmers}){
				#warn substr($sKmer,$i,1),"\n";
				if (substr( $sKmer, $i , 1 ) eq $b){
					$rPFM->[$i]->{$b}+=$raZscores->[$n];
					#warn "".$rPFM->[$i]->{$b}."\n";
					$anTotals[$i]+=$raZscores->[$n];
				}elsif( substr( $sKmer, $i , 1 ) eq '-') {
					$rPFM->[$i]->{$b}+=($avgZ/4);
					$anTotals[$i]+=$avgZ;
				}
				$n++;
			}
		}
	}
	
	# convert to probabilities
	
	#warn join(";",@anTotals),"\n";
	
	for(my $i=0; $i<$nWidth; $i++){ #position in width
		foreach my $b (@aAlphabet){
			if($bSmallSampleCor){
				$rPFM->[$i]->{$b} = ( $rPFM->[$i]->{$b} + sqrt($anTotals[$i])*0.25 ) / ($anTotals[$i] + sqrt($anTotals[$i]));
			} elsif (!$bGetCounts) {
				$rPFM->[$i]->{$b} = $rPFM->[$i]->{$b} / $anTotals[$i];
			}
		}
	}
	
	return $rPFM;
	#print "PFM!\n";
	#pretty_print_pssm($rPFM);

}


__DATA__

pfm_from_aligned.pl [FILE | < FILE]

Reads in a file (or piped input) containing a list of aligned k-mers and
spits out a PFM.

Note 1: gaps ('-') characters give a weight of 0.25 to all bases. This is
primarily so that logos made from this PFM don't have high-stringency bases
on the edges. TODO: make a flag to adjust this behaviour.

Note 2: RNA-based alphabet.

   -s:  apply small sample correction TODO: show math here
   		Note: incompatible with the -c option
   -c:  use frequency counts instead of fractions (ie dont divide by the number of sites)
   -g:	don't set '-' characters to 0.25 for all bases (allows for trimming later)

