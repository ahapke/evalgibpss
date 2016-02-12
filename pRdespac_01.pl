#!/usr/bin/perl -w
#pRdespac version 01 Copyright 2015 Andreas Hapke
#This program reads the file *.alleles produced by pyRAD version 3.0.1
#removes the spacer and surrounding gaps from each sequence in the file
#and prints the new sequences to an outfile.
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#Contact:
#Andreas Hapke 
#Johannes Gutenberg University Mainz
#Institute of Anthropology
#Anselm-Franz-von-Benzel-Weg 7
#D-55099 Mainz, Germany
#E-Mail: ahapke2@gmail.com

use strict;
use warnings;

#######################
#Declare and initialize
#######################
my $infn = '';
my $inline = '';
my $inlineNo = 0;
my @inloc = ();#all sequence containing lines of one locus
my $ngaps1 = 0;#number of gaps before spacer
my $ngaps2 = 0;#number of gaps after spacer
my @ngapsleft = ();#numbers of gaps before spacer in one locus
my @ngapsright = ();#numbers of gaps after spacer in one locus
my $maxgaps1 = 0;#greatest number of gaps before spacer in one locus
my $maxgaps2 = 0;#greatest number of gaps after spacer in one locus
my $spacerstart = 0;#start position of spacer
my %spacerstarts = ();#collects start positions of spacer in one loc for check
my $nspacerstarts = 0;#number of spacer start positions in a locus, must be 1
my $splicestart = 0;#start position of alignment part to remove
my $splicelen = 0;# length of alignment part to remove
my $outfn = '';
my $outline = '';

#open infile
$infn = $ARGV[0];
unless(open(INFILE, "<", $infn)) {
	print "Cannot open $infn, exiting..\n";
}
$outfn = 'nospac_' . $infn;
unless(open(OUTFILE, ">", $outfn)) {
	print "Cannot open $outfn, exiting..\n";
	exit;
}

#read infile locus-wise and treat
while ($inline = <INFILE>) {
	chomp $inline;
	++$inlineNo;
	if ($inline =~ /(-*)(nnnn)(-*)/g) {
		$ngaps1 = length $1;
		$ngaps2 = length $3;
		$spacerstart = pos($inline) - $ngaps2 - 4;
		push @inloc, $inline;
		push @ngapsleft, $ngaps1;
		push @ngapsright, $ngaps2;
		++$spacerstarts{$spacerstart};		
	}
	elsif ($inline =~ m{^//}) {#reached last line of locus
		if (keys %spacerstarts > 1) {#not all spacers at same start pos? Exit!
			print "Locus ending at line $inlineNo: spacers not aligned!, exiting..\n";
			exit;
		}
		@ngapsleft = reverse sort {$a <=> $b} @ngapsleft;
		@ngapsright = reverse sort {$a <=> $b} @ngapsright;
		$maxgaps1 = $ngapsleft[0];
		$maxgaps2 = $ngapsright[0];
		$spacerstart = (sort {$a <=> $b} keys %spacerstarts)[0];
		$splicestart = $spacerstart - $maxgaps1;
		$splicelen = $maxgaps1 + 4 + $maxgaps2;
		#remove middle part from each line and print to outfile
		for $outline (@inloc) {
			substr($outline, $splicestart, $splicelen) = '';
			print OUTFILE "$outline\n";
		}
		#Don't forget the last line:
		$outline = $inline;
		substr($outline, $splicestart, $splicelen) = '';
		print OUTFILE "$outline\n";
		#set variables back and next loc
		@inloc = ();
		@ngapsleft = ();
		@ngapsright = ();
		%spacerstarts = ();
	}
	else {
		print "Line $inlineNo: I don't understand this line, exiting..\n";
		exit;
	}
}
close INFILE;
close OUTFILE;
