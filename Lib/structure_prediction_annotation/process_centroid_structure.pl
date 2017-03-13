#!/usr/bin/perl

use strict;
use warnings;

my $id = '';
my $centroid_struct = '';

while(<>){
	if (/^>/){
		my $newid = $_;
		print "$id$centroid_struct\n";
		$id = $newid;
	} elsif (/^(.+)\s+{.+}$/) {
		$centroid_struct = $1;
	}
}
print "$id$centroid_struct\n";