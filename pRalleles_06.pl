#!/usr/bin/perl -w
#pRalleles version 06 Copyright 2015 Andreas Hapke
#This program analyzes locus coverage and performs pairwise comparisons
#of multi-locus genotypes based on output of pyRAD version 3.0.1.
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

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#######################
#Declare and initialize
#######################

#{
my $pairfn = '';#name of file with pairs of indIDs to compare, separated by TAB, one pair per line
my $infn = '';#name of pyRAD *.alleles file
my $outfn = '';#name of an outfile to be produced by a sub
my %pairs = ();# {pairNo}=(ind1,ind2)
my %nonidpairs = ();# {pairNo} = (ind1, ind2) all pairs of inds not defined in user file
my $pairNo = 0;
my %ingt = ();# {loc}{ind}{allNo} = seq
my %indloc = ();# {ind} = number of loci
my $loc = 0;
my $ind = '';
my $ind1 = '';
my $ind2 = '';
my $allNo = '';
my $seq = '';
my $tempstring1 = '';
my $tempstring2 = '';
my @temparr1 = ();
my %temphash1 = ();
my $i = 0;
my $j = 0;
my $k = 0;
#}

##############
#Read in pairs
##############

#{
$pairfn = shift @ARGV;
unless(open(INFILE, "<", $pairfn)) {
 print "Cannot open $pairfn, exiting..\n";
 exit;
}
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split (/\t/, $tempstring1);
	@{$pairs{$pairNo}} = @temparr1;
	++$pairNo;
}
close INFILE;
#}

########################
#Read in *.alleles files
#Populate %ingt
#Populate %indloc
########################

#{
while (@ARGV) {
	$infn = shift @ARGV;
	unless(open(INFILE, "<", $infn)) {
	 print "Cannot open $infn, exiting..\n";
	 exit;
	}
	while ($tempstring1 = <INFILE>) {
		chomp $tempstring1;
		if ($tempstring1 =~ /^>/) {
			$tempstring1 =~ s/^>//;
			@temparr1 = split (/\s+/, $tempstring1);
			$seq = $temparr1[1];
			$tempstring2 = $temparr1[0];
			@temparr1 = split(/_/,$tempstring2);
			$ind = $temparr1[0];
			$ind =~ s/.assembled$//;
			$allNo = $temparr1[1];
			$ingt{$loc}{$ind}{$allNo} = $seq;
			$indloc{$ind} += 0.5;
		}
		elsif ($tempstring1 =~ m{^//}) {
			++$loc;
		}
	}
	close INFILE;
}
#}

###########################
#Produce outfile indloc.txt
###########################

#{
unless(open(OUTFILE, ">", "indloc.txt")) {
 print "Cannot open outfile indloc.txt, exiting..\n";
 exit;
}
print OUTFILE "ind\tnloci\n";
for $ind (sort keys %indloc) {
	print OUTFILE "$ind\t$indloc{$ind}\n";
}
close OUTFILE;
#}

#####################
#Populate %nonidpairs
#####################

#{
#populate a hash with all pairs defined in user file as keys, in both orders per pair
for $pairNo (keys %pairs) {
	$ind1 = $pairs{$pairNo}[0];
	$ind2 = $pairs{$pairNo}[1];
	$tempstring1 = $ind1 . '_' . $ind2;
	$temphash1{$tempstring1} = 1;
	$tempstring1 = $ind2 . '_' . $ind1;
	$temphash1{$tempstring1} = 1;
}
#identify all pairs not defined in user files and populate %nonidpairs
@temparr1 = keys %indloc;
@temparr1 = sort @temparr1;
$k = @temparr1;
$pairNo = 0;
for ($i = 0; $i < $k-1; ++$i) {
	for ($j = $i+1; $j < $k; ++$j) {
		$ind1 = $temparr1[$i];
		$ind2 = $temparr1[$j];
		$tempstring1 = $ind1 . '_' . $ind2;
		$tempstring2 = $ind2 . '_' . $ind1;
		unless((defined $temphash1{$tempstring1}) or (defined $temphash1{$tempstring2})) {
			@{$nonidpairs{$pairNo}} = ($ind1,$ind2);
			++$pairNo;
		}
	}
}
#}

#############
#Analyze data
#############

#{
#Call sub coverage to analyze number of loci and individuals
coverage(\%ingt);
#Compare pairs of individuals defined in user file
$outfn = 'idpair_comp.txt';
indpairs(\%pairs,\%ingt,$outfn);
#Compare pairs of individuals NOT defined in user file
$outfn = 'nonidpair_comp.txt';
indpairs(\%nonidpairs,\%ingt,$outfn);
#Count genotypes with indel variation
countindelgt(\%ingt);
#Count loci with indel variation
countindelloc(\%ingt);
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";

exit;


############
#Subroutines
############

#Definition of subroutine coverage
#Expects:
#ref to %ingt {loc}{ind}{allNo} = seq
#Counts loci and number of individuals for each locus
#Prints counts to an outfile coverage.txt
sub coverage {
	#declare and initialize: _r means a reference
	my($ingt_r) = @_;
	my $loc = 0;
	my $ind = '';
	my $nind = 0;
	my $nindmax = 0;#maximum number of individuals in a locus
	my $allNo = 0;
	my %coverage = ();# {nInd} = nLoc
	my %nloccum = ();# {nInd} = nLoc cumulative: loci with AT LEAST 1,2,3,4... inds
	my %nloccumperc = ();# nloccum as percent of all loci
	my $nloc = 0;
	my $i = 0;
	my $j = 0;
	
	#determine number of loci with exactly n inds
	for $loc (keys %{$ingt_r}) {
		for $ind (keys %{$$ingt_r{$loc}}) {
			++$nind;
		}
		++$coverage{$nind};
		$nind = 0;
	}
	#set nonexisting counts of individuals to 0
	$nindmax = (reverse sort {$a <=> $b} keys %coverage)[0];
	for ($i = 0; $i <= $nindmax; ++$i) {
		unless(defined $coverage{$i}) {
			$coverage{$i} = 0;
		}
	}
	
	#determine cumulative locus counts
	for ($i = 1; $i <= $nindmax; ++$i) {
		for ($j = $i; $j <= $nindmax; ++ $j) {
			if (defined $coverage{$j}) {
				$nloc += $coverage{$j};
			}
		}
		$nloccum{$i} = $nloc;
		$nloc = 0;
	}
	#determine cumulative percentages
	for ($i = 1; $i <= $nindmax; ++$i) {
		$nloccumperc{$i} = $nloccum{$i} / $nloccum{1} * 100;
	}
	
	#produce outfile coverage.txt
	unless(open(OUTFILE, ">", "coverage.txt")) {
		print "Cannot open outfile coverage.txt, exiting..\n";
		exit;
	}
	print OUTFILE "nInd\tnloc\tnloccum\tnloccumperc\n";
	for ($i = 1; $i <= $nindmax; ++$i) {
		print OUTFILE "$i\t$coverage{$i}\t$nloccum{$i}\t$nloccumperc{$i}\n";
	}
	close OUTFILE;
}

#Definition of subroutine indpairs
#Expects:
#ref to %pairs: {pairNo} = (ind1,ind2)
#ref to %ingt: {loc}{ind}{allNo} = seq
#$outfn: name of outfile to produce
#
#Performs pairwise comparisons of individuals
#nloc1: number of loci in ind1
#nloc2: number of loci in ind2
#nlocboth: number of loci in typed in both individuals
#nloctot: total number of loci in both individuals
#weakpboth: nlocboth/loci in stronger ind *100
#bothptot: nlocboth / nloctotal * 100
#nlocid: number of genotypes identical in both
#percid: percent of loci identical in both (relative to nlocboth)
#identical: N counts as identical with -ACGT
#min, max, and average of these
sub indpairs {
	#declare and initialize: _r means a reference
	my($pairs_r, $ingt_r, $outfn) = @_;
	my $pairNo = 0;#number of a pair
	my $ind1 = '';# ID of ind1 of a pair
	my $ind2 = '';# ID of ind2 of a pair
	my $ind1seq1 = '';# seq1 of ind1 at a locus
	my $ind1seq2 = '';
	my $ind2seq1 = '';
	my $ind2seq2 = '';
	my $dist = 0;
	my $nlocind1 = 0;# number of loci in ind1
	my $nlocind2 = 0;# number of loci in ind2
	my $nlocboth = 0;# number of loci in both
	my $nloctot = 0;# total number of loci in both inds
	my $nlocweak = 0;# number of loci in weaker ind
	my $weakpboth = 0;# nlocboth / nloc weaker ind *100
	my $bothptot = 0;# nlocboth / nloctot * 100
	my $nlocid = 0;# number of loci with identical genotypes in both
	my $percid = 0;# percent of loci with identical genotypes in both
	my $loc = 0;# ID of a locus
	my %pairres = ();#results for each pair:
					# {pairNo} = (nloc1,nloc2,nlocboth,nloctot,weakpstrong,bothptot,nlocid,perclocid);
	my $nres = 8;#number of result values per pair
	my $npairs = 0;#number of ind pairs
	my $total = 0;#sum over result values
	my $min = 0;#min
	my $max = 0;#max
	my $avg = 0;#avg
	my @avgres = ();#averages over pairs as in %pairres
	my @minres = ();#min over pairs as in %pairres
	my @maxres = ();#max over pairs as in %pairres
	my @temparr1 = ();
	my $tempstring1 = '';
	my $i = 0;
	
	#open outfile and print header line
	unless(open(OUTFILE, ">", $outfn)) {
		print "Cannot open $outfn, exiting..\n";
		exit;
	}
	print OUTFILE "ind1\tind2\tnloc1\tnloc2\tnlocboth\tnloctot\tpwnest\tpwover\tnlocid\tperclocid\n";
	
	for $pairNo (sort {$a <=> $b} keys %{$pairs_r}) {#loop through indpairs
		$ind1 = $$pairs_r{$pairNo}[0];
		$ind2 = $$pairs_r{$pairNo}[1];
		for $loc (keys %{$ingt_r}) {
			if (defined $$ingt_r{$loc}{$ind1}) {
				++$nlocind1;
			}
			if (defined $$ingt_r{$loc}{$ind2}) {
				++$nlocind2;
			}
			if ((defined $$ingt_r{$loc}{$ind1}) and (defined $$ingt_r{$loc}{$ind2})) {
				#locus typed in both individuals
				++$nlocboth;
				$ind1seq1 = $$ingt_r{$loc}{$ind1}{0};
				$ind1seq2 = $$ingt_r{$loc}{$ind1}{1};
				$ind2seq1 = $$ingt_r{$loc}{$ind2}{0};
				$ind2seq2 = $$ingt_r{$loc}{$ind2}{1};
				#compare the genotypes
				$dist = distpwdN($ind1seq1,$ind2seq1);
				if ($dist == 0) {
					$dist = distpwdN($ind1seq2,$ind2seq2);
					if ($dist == 0) {
						++$nlocid;
					}
				} else {
					$dist = distpwdN($ind1seq2, $ind2seq1);
					if ($dist == 0) {
						$dist = distpwdN($ind1seq1, $ind2seq2);
						if ($dist == 0) {
							++$nlocid;
						}
					}
				}
			}
		}
		#calculate results for this pair and store in %pairres
		$nloctot = $nlocind1 + $nlocind2 - $nlocboth;
		@temparr1 = ($nlocind1,$nlocind2);
		@temparr1 = sort {$a <=> $b} @temparr1;
		$nlocweak = $temparr1[0];
		$weakpboth = $nlocboth / $nlocweak * 100;
		$bothptot = $nlocboth / $nloctot * 100;
		$percid = $nlocid / $nlocboth * 100;
		@{$pairres{$pairNo}} = ($nlocind1,$nlocind2,$nlocboth,$nloctot,$weakpboth,$bothptot,$nlocid,$percid);
		#set variables back and next ind
		$nlocind1 = 0;
		$nlocind2 = 0;
		$nlocboth = 0;
		$nloctot = 0;
		$nlocweak = 0;
		$weakpboth = 0;
		$bothptot = 0;
		$nlocid = 0;
		$percid = 0;
	}
	#calculate average, min, max
	$npairs = keys %pairres;
	@temparr1 = ();
	for ($i = 0; $i < $nres; ++$i) {
		for $pairNo (sort {$a <=> $b} keys %pairres) {
			$total += $pairres{$pairNo}[$i];
			push @temparr1, $pairres{$pairNo}[$i];
		}
		@temparr1 = sort {$a <=> $b} @temparr1;
		$avg = $total / $npairs;
		$min = $temparr1[0];
		$max = pop @temparr1;
		push @avgres, $avg;
		push @minres, $min;
		push @maxres, $max;
		#set variables back and next i
		$total = 0;
		@temparr1 = ();
	}
	#print results to outfile
	$tempstring1 = join("\t",@avgres);
	print OUTFILE "avg\tavg\t$tempstring1\n";
	$tempstring1 = join("\t",@minres);
	print OUTFILE "min\tmin\t$tempstring1\n";
	$tempstring1 = join("\t",@maxres);
	print OUTFILE "max\tmax\t$tempstring1\n";
	for $pairNo (sort {$a <=> $b} keys %pairres) {
		$ind1 = $$pairs_r{$pairNo}[0];
		$ind2 = $$pairs_r{$pairNo}[1];
		$tempstring1 = join("\t",@{$pairres{$pairNo}});
		print OUTFILE "$ind1\t$ind2\t$tempstring1\n";
	}	
	close OUTFILE;
}

#definition of subroutine distpwdN version 05
#calculates number of differences between 2 DNA sequences of equal length
#with pairwise deletion of positions containing N
#knows chars: : -ACGTN
#expects 2 scalars
#$seq1, $seq2
#returns $dist
sub distpwdN {
	#declare and initialize
	my ($seq1,$seq2) = @_;
	my $seql = 0;# sequence length
	my $dist = 0;#the result: number of differences
	my $Nstring = '';#a string of N of length seqlength
	my $NAseq1 = '';#seq1 with CGT replaced by AAA
	my $NAseq2 = '';#seq2 with CGT replaced by AAA
	my $changepos = '';#string that contains defined values at positions
						#where only one of the 2 seqs has N, other pos: \0
	my$pos = 0;#a position;
	
	$seql = length $seq1;
	$Nstring = "N" x $seql;#create string of seqlength x N
	#create the 2 NAseqs
	($NAseq1 = $seq1) =~ tr/-CGT/AAAA/;
	($NAseq2 = $seq2) =~ tr/-CGT/AAAA/;
	#determine N-positions in seq1 and seq2 as \0 others as defined
	#match these against each other to get $changepos
	#:pos with N in only 1 seq: defined, others \0
	$changepos = (($NAseq1 ^ $Nstring) ^ ($NAseq2 ^ $Nstring));
	#loop through $changepos, for each changepos, replace corresponding
	#character in seq1 and seq2 with N
	while ($changepos =~ /[^\0]/g) {
		$pos = sprintf "%d", pos($changepos)-1;
		substr($seq1,$pos,1) = 'N';
		substr($seq2,$pos,1) = 'N';
	}
	#calculate the distance
	$dist = $seql - (($seq1 ^ $seq2) =~ tr/\0//);
	return $dist;	
}
#definition of sub countindelgt
#Expects:
#ref to %ingt {loc}{ind}{allNo} = seq
#Analyzes all individual genotypes
#Counts genotypes with indel variation between the two alleles
#Produces outfile indelgtcount.txt
sub countindelgt {
	my ($ingt_r) = @_;
	my $loc = 0;
	my $ind = 0;
	my $seq1 = '';#seq1 of a genotype
	my $seq2 = '';#seq2 of a genotype
	my $nindelgt = 0;#number of gentoypes with indel variation
	my $dist = 0;
	
	#open outfile
	unless(open(OUTFILE, ">", "indelgtcount.txt")) {
		print "Cannot open outfile indelgtcount.txt, exiting..\n";
		exit;
	}
	
	for $loc (keys %{$ingt_r}) {
		for $ind (keys %{$$ingt_r{$loc}}) {
			$seq1 = $$ingt_r{$loc}{$ind}{0};
			$seq2 = $$ingt_r{$loc}{$ind}{1};
			#if the 2 seqs differ and at least one contains at least one gap
			if (($seq1 ne $seq2) and (($seq1 =~ /-/) or $seq2 =~ /-/)) {
				#change all characters except "-" to A
				$seq1 =~ tr/NCGT/AAAA/;
				$seq2 =~ tr/NCGT/AAAA/;
				#calculate distance with distpwdN
				$dist = distpwdN($seq1,$seq2);
				if ($dist > 0) {
					++$nindelgt;
				}
			}
		}
	}
	print OUTFILE "number of genotypes with indel variation:\t$nindelgt\n";
	close OUTFILE;
}

#definition of sub countindelloc
#Expects:
#ref to %ingt {loc}{ind}{allNo} = seq
#Counts loci with indel variation and total number of indels
#Counts every indel position as one indel
#Distinguishes between two kinds of indels:
#Positions that contain only - and N: N-indel
#Positions that contain ACGT and - indel
#Produces outfile indel_count.txt
sub countindelloc {
	#declare and initialize: _r means a reference
	my($ingt_r) = @_;
	my $loc = 0;
	my $ind = '';
	my $seq = '';
	my %seqvars = ();# {seq} = count collapses sequences into nonidentical variants
	my @ali = ();# 2d array, contains all nonidentical sequences of one loc, d1: rows, d2 cols
	my $hasgap = 0;# 1 if a locus contains at least one gap
	my $col = 0;#column in a locus alignment
	my $row = 0;#row in a locus alignment
	my %charcount = ();# {char} = count counts characters in a column of a locus alignment
	my $nNindelthisloc = 0;# number of indel positions in a locus that contain only - and N
	my $nindelthisloc = 0;#number of indel positions in a locus that contain - and at least one of ACGT
	my $nNindel = 0;# total number of N- indel positions
	my $nindel = 0;#total number of indel positions
	my $nNindelloc = 0;#total number of loci with at least one N- indel but no -/ACGT indel
	my $nindelloc = 0;#total number of loci with at least one indel
	my $percNindelloc = 0;#percent of loci with at least one N- indel but no -/ACGT indel
	my $percindelloc = 0;#percent of loci with at least one indel
	my $nloc = 0;#total number of loci
	my @temparr1 = ();
	
	#open outfile and print header line
	unless(open(OUTFILE, ">", "indel_count.txt")) {
		print "Cannot open outfile indel_count.txt, exiting..\n";
		exit;
	}
	print OUTFILE "nlociwithindel\tnlociwithonlyN-indel\tperclociwithindel\tperclociwithonlyN-indel\tnindel\tnN-indel\n";
	
	for $loc (keys %{$ingt_r}) {
		++$nloc;
		for $ind (keys %{$$ingt_r{$loc}}) {
			$seq = $$ingt_r{$loc}{$ind}{0};
			++$seqvars{$seq};#collapse sequences into nonidentical ones as keys in %seqvars
			$seq = $$ingt_r{$loc}{$ind}{1};
			++$seqvars{$seq};#collapse sequences into nonidentical ones as keys in %seqvars
		}
		#Check if locus contains indel variation
		for $seq (keys %seqvars) {
			if ($seq =~ /-/) {
				$hasgap = 1;
				last;
			}
		}
		if ($hasgap == 1) {#if locus contains indel variation
			#build 2d array @ali: d1 rows, d2 cols
			for $seq (keys %seqvars) {
				@temparr1 = split('',$seq);
				push @ali, [@temparr1];
			}
			#find indel positions in @ali
			for ($col = 0; $col < @{$ali[0]}; ++$col) {
				for ($row = 0; $row < @ali; ++$row) {
					++$charcount{${$ali[$row]}[$col]};
				}
				if ((defined $charcount{"-"}) and (keys %charcount >= 2)) {#if column contains gap and at least one other character
					if ((defined $charcount{"N"}) and (keys %charcount == 2)) {#if columns contains only gap and N
						++$nindelthisloc;#count an indel
						++$nNindelthisloc;#count an N-indel
					} else {#if column contains gap and at least one of ACGT
						++$nindelthisloc;#count a "good" indel
					}
				}
				%charcount = ();#set back and next col
			}
			#count indels
			$nindel += $nindelthisloc;
			$nNindel += $nNindelthisloc;
			#count loci with indels
			if ($nindelthisloc > 0) {#if locus has any indel
				++$nindelloc;#count
				if ($nindelthisloc == $nNindelthisloc) {#if locus has only N- indels, count 
					++$nNindelloc;#count
				}
			}
		}
		#set variables back and next loc
		%seqvars = ();
		$hasgap = 0;
		@ali = ();
		$nNindelthisloc = 0;
		$nindelthisloc = 0;
	}
	#calculate percentages
	$percindelloc = $nindelloc / $nloc * 100;
	$percNindelloc = $nNindelloc / $nloc *100;
	#print results to outfile
	print OUTFILE "$nindelloc\t$nNindelloc\t$percindelloc\t$percNindelloc\t$nindel\t$nNindel\n";
	close OUTFILE;
}