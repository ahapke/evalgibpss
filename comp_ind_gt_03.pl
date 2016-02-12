#!/usr/bin/perl -w
#comp_ind_gt version 03 Copyright 2015 Andreas Hapke
#This program performs pairwise comparisons of multi-locus genotypes
#based on output of GIbPSs.
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
#declare and initialize
#######################

#{
my $pairsfilename = '';#file with IDs of paired inds, one pair per line
my $genotypesfilename = 'export/genotypes.txt';#genotypes file out of export_tables
my $sumoutname = 'export/comp_ind_gt_sum.txt';#outfile with summary counts
my %ind_col = ();# {ind}=column index for genotypes array in %genotypes
my %genotypes = ();# {poplocID}=@ of genotypes as in file genotypes.txt
my $nloc_total = 0;#total number of loci in dataset
my $ind = '';#ind ID
my $ind1 = '';#ind1 of a pair
my $ind2 = '';#ind2 of a pair
my $ind1col = 0;#index for genotype
my $ind2col = 0;#index for genotype
my $all_1 = '';#allele of ind1
my $all_2 = '';#allele of ind2
my $nall_1 = 0;#number of alleles in ind1
my $nall_2 = 0;#number of alleles in ind2
my %alleles1 = ();#all alleles of ind1
my %alleles2 = ();#all alleles of ind2
my $gt1 = '';#genotype1
my $gt2 = '';#genotype2
my $n_mis = 0;#number of mismatches
my $poplocID = 0;
my %pairs = ();# {pair_No} = (ind1,ind2)
my $pair_No = 0;
my %counts = ();#counts results
my $tempstring1 = '';
my @temparr1 = ();
my $i = 0;
my $val = 0;
#}

#############
#Preparations
#############

#{
#read in genotypes
unless(open(INFILE,$genotypesfilename)) {
	print "Cannot open $genotypesfilename Exiting..\n";
	exit;
}
#get individual column indices out of headerline
$tempstring1 = <INFILE>;
chomp $tempstring1;
@temparr1 = split(/\t/,$tempstring1);
shift @temparr1;#remove first col
$i = 0;
for $ind (@temparr1) {
	$ind_col{$ind} = $i;
	++$i;
}
#read in genotypes
while($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	$poplocID = shift @temparr1;
	@{$genotypes{$poplocID}} = @temparr1;
	++$nloc_total;
}
close INFILE;
#ask user for pairs file and open
print "I need a file with pairs of ind IDs:\n",
"two IDs per line, separated by TAB.\n",
"I will preserve the order of inds in each line.\n",
"Please enter filename: ";
$pairsfilename = <STDIN>;
chomp $pairsfilename;
unless(open(INFILE,$pairsfilename)) {
	print "Cannot open $pairsfilename Exiting..\n";
	exit;
}
#read in pairs file
while($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	$tempstring1 =~ s/^\s*//;#remove leading whitespace
	$tempstring1 =~ s/\s*$//;#remove trailing whitespace
	@temparr1 = split(/\t/,$tempstring1);
	unless(@temparr1 == 2) {
		print "I do not understand this line:\n$tempstring1\nExiting..\n";
		exit;
	}
	$ind1 = $temparr1[0];
	$ind2 = $temparr1[1];
	unless(defined $ind_col{$ind1}) {
		print "I do not know individual $ind1. Exiting..\n";
	}
	unless(defined $ind_col{$ind2}) {
		print "I do not know individual $ind2. Exiting..\n";
	}	
	@{$pairs{$pair_No}} = ($ind1,$ind2);
	++$pair_No;
}
close INFILE;

#initialize %counts
for $pair_No (keys %pairs) {
	$counts{'nloc_ind1'}{$pair_No} = 0;
	$counts{'nloc_ind2'}{$pair_No} = 0;
	$counts{'nloc_both'}{$pair_No} = 0;
	$counts{'id'}{$pair_No} = 0;
	$counts{'drop'}{$pair_No} = 0;
	$counts{'mis'}{$pair_No} = 0;
	$counts{'id_homhom'}{$pair_No} = 0;
	$counts{'id_hethet'}{$pair_No} = 0;
	$counts{'id_other'}{$pair_No} = 0;
	$counts{'drop_homhet'}{$pair_No} = 0;
	$counts{'drop_hethom'}{$pair_No} = 0;
	$counts{'drop_other'}{$pair_No} = 0;
	$counts{'mis_homhom'}{$pair_No} = 0;
	$counts{'mis_hethet'}{$pair_No} = 0;
	$counts{'mis_homhet'}{$pair_No} = 0;
	$counts{'mis_hethom'}{$pair_No} = 0;
	$counts{'mis_other'}{$pair_No} = 0;
	$counts{'NN_ind1'}{$pair_No} = 0;
	$counts{'NN_ind2'}{$pair_No} = 0;
	$counts{'NN_both'}{$pair_No} = 0;
}
#create outfile
unless(open(SUMOUT, ">$sumoutname")) {
	print "Cannot open $sumoutname. Exiting\n";
	exit;
}
#print headerline
print SUMOUT "category";
for $pair_No (sort {$a <=> $b} keys %pairs) {
	print SUMOUT "\t$pairs{$pair_No}[0]","_","$pairs{$pair_No}[1]";
}
print SUMOUT "\n";	
#}

##################
#Compare genotypes
##################

#{
for $pair_No (keys %pairs) {
	$ind1 = $pairs{$pair_No}[0];
	$ind2 = $pairs{$pair_No}[1];
	$ind1col = $ind_col{$ind1};
	$ind2col = $ind_col{$ind2};
	for $poplocID (keys %genotypes) {
		$gt1 = $genotypes{$poplocID}[$ind1col];
		$gt2 = $genotypes{$poplocID}[$ind2col];
		#Genotype missing in both?
		if (($gt1 eq '-999') and ($gt2 eq '-999')) {
			++$counts{'NN_both'}{$pair_No};
			++$counts{'NN_ind1'}{$pair_No};
			++$counts{'NN_ind2'}{$pair_No};
		}
		#Genotype missing in ind1?
		elsif ($gt1 eq '-999') {
			++$counts{'NN_ind1'}{$pair_No};
			++$counts{'nloc_ind2'}{$pair_No};
		}
		#Genotype missing in ind2?
		elsif ($gt2 eq '-999') {
			++$counts{'NN_ind2'}{$pair_No};
			++$counts{'nloc_ind1'}{$pair_No};
		}
		else {#Compare genotypes
			++$counts{'nloc_ind1'}{$pair_No};
			++$counts{'nloc_ind2'}{$pair_No};
			++$counts{'nloc_both'}{$pair_No};
			#get alleles of ind1
			@temparr1 = split(/\//,$gt1);
			$nall_1 = @temparr1;
			for $all_1 (@temparr1) {
				$alleles1{$all_1} = 1;
			}
			#get alleles of ind2
			@temparr1 = split(/\//,$gt2);
			$nall_2 = @temparr1;
			for $all_2 (@temparr1) {
				$alleles2{$all_2} = 1;
			}
			if ($nall_1 == $nall_2) {#Number of alleles equal
				$n_mis = 0;
				for $all_1 (keys %alleles1) {#count mismatches
					unless (defined $alleles2{$all_1}) {
						++$n_mis;
					}
				}
				if ($n_mis == 0) {#all alleles identical
					++$counts{'id'}{$pair_No};
					if ($nall_1 == 1) {#id_homhom
						++$counts{'id_homhom'}{$pair_No};
					}
					elsif ($nall_1 == 2) {#id_hethet
						++$counts{'id_hethet'}{$pair_No};
					}
					else {#id more alleles
						++$counts{'id_other'}{$pair_No};
					}
				} else {#mismatch
					++$counts{'mis'}{$pair_No};
					if ($nall_1 == 1) {#mis_homhom
						++$counts{'mis_homhom'}{$pair_No};
					}
					elsif ($nall_1 == 2) {#mis_hethet
						++$counts{'mis_hethet'}{$pair_No};
					}
					else {#mis more alleles
						++$counts{'mis_other'}{$pair_No};
					}					
				}
			} else {#Number of alleles different
				if ($nall_1 < $nall_2) {#ind 1 has less alleles
					$n_mis = 0;
					for $all_1 (keys %alleles1) {#count mismatches
						unless (defined $alleles2{$all_1}) {
							++$n_mis;
						}
					}
					if ($n_mis == 0) {#all alleles of ind1 match ind2
						++$counts{'drop'}{$pair_No};
						if (($nall_1 == 1) and ($nall_2 == 2)) {#drop homhet
							++$counts{'drop_homhet'}{$pair_No};
						} else {#drop other
							++$counts{'drop_other'}{$pair_No};
						}					
					} else {#mismatch
						++$counts{'mis'}{$pair_No};
						if (($nall_1 == 1) and ($nall_2 == 2)) {#mis homhet
							++$counts{'mis_homhet'}{$pair_No};
						} else {#mis other
							++$counts{'mis_other'}{$pair_No};
						}						
					}
				}
				if ($nall_1 > $nall_2) {#ind 2 has less alleles
					$n_mis = 0;
					for $all_2 (keys %alleles2) {#count mismatches
						unless (defined $alleles1{$all_2}) {
							++$n_mis;
						}
					}
					if ($n_mis == 0) {#all alleles of ind2 match ind1
						++$counts{'drop'}{$pair_No};
						if (($nall_1 == 2) and ($nall_2 == 1)) {#drop hethom
							++$counts{'drop_hethom'}{$pair_No};
						} else {#drop other
							++$counts{'drop_other'}{$pair_No};
						}					
					} else {#mismatch
						++$counts{'mis'}{$pair_No};
						if (($nall_1 == 2) and ($nall_2 == 1)) {#mis hethom
							++$counts{'mis_hethom'}{$pair_No};
						} else {#mis other
							++$counts{'mis_other'}{$pair_No};
						}						
					}
				}			
			}
		}
		%alleles1 = ();#set back
		%alleles2 = ();#set back and next poploc
	}
}
#}

##############################
#Print data to summary outfile
##############################

#{
print SUMOUT "nloc_total";
for $pair_No (sort {$a <=> $b} keys %{$counts{'nloc_ind1'}}) {
	print SUMOUT "\t$nloc_total";
}
print SUMOUT "\n";
print SUMOUT "nloc_ind1";
for $pair_No (sort {$a <=> $b} keys %{$counts{'nloc_ind1'}}) {
	print SUMOUT "\t$counts{'nloc_ind1'}{$pair_No}";
}
print SUMOUT "\n";
print SUMOUT "nloc_ind2";
for $pair_No (sort {$a <=> $b} keys %{$counts{'nloc_ind2'}}) {
	print SUMOUT "\t$counts{'nloc_ind2'}{$pair_No}";
}
print SUMOUT "\n";
print SUMOUT "nloc_both";
for $pair_No (sort {$a <=> $b} keys %{$counts{'nloc_both'}}) {
	print SUMOUT "\t$counts{'nloc_both'}{$pair_No}";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: identical";
for $pair_No (sort {$a <=> $b} keys %{$counts{'id'}}) {
	$val = $counts{'id'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: dropout";
for $pair_No (sort {$a <=> $b} keys %{$counts{'drop'}}) {
	$val = $counts{'drop'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: mismatch";
for $pair_No (sort {$a <=> $b} keys %{$counts{'mis'}}) {
	$val = $counts{'mis'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: id_homhom";
for $pair_No (sort {$a <=> $b} keys %{$counts{'id_homhom'}}) {
	$val = $counts{'id_homhom'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: id_hethet";
for $pair_No (sort {$a <=> $b} keys %{$counts{'id_hethet'}}) {
	$val = $counts{'id_hethet'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: id_other";
for $pair_No (sort {$a <=> $b} keys %{$counts{'id_other'}}) {
	$val = $counts{'id_other'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: drop_homhet";
for $pair_No (sort {$a <=> $b} keys %{$counts{'drop_homhet'}}) {
	$val = $counts{'drop_homhet'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: drop_hethom";
for $pair_No (sort {$a <=> $b} keys %{$counts{'drop_hethom'}}) {
	$val = $counts{'drop_hethom'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: drop_other";
for $pair_No (sort {$a <=> $b} keys %{$counts{'drop_other'}}) {
	$val = $counts{'drop_other'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: mis_homhom";
for $pair_No (sort {$a <=> $b} keys %{$counts{'mis_homhom'}}) {
	$val = $counts{'mis_homhom'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: mis_hethet";
for $pair_No (sort {$a <=> $b} keys %{$counts{'mis_hethet'}}) {
	$val = $counts{'mis_hethet'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: mis_homhet";
for $pair_No (sort {$a <=> $b} keys %{$counts{'mis_homhet'}}) {
	$val = $counts{'mis_homhet'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: mis_hethom";
for $pair_No (sort {$a <=> $b} keys %{$counts{'mis_hethom'}}) {
	$val = $counts{'mis_hethom'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of nloc_both: mis_other";
for $pair_No (sort {$a <=> $b} keys %{$counts{'mis_other'}}) {
	$val = $counts{'mis_other'}{$pair_No} / $counts{'nloc_both'}{$pair_No} * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of all loci: missing ind 1";
for $pair_No (sort {$a <=> $b} keys %{$counts{'NN_ind1'}}) {
	$val = $counts{'NN_ind1'}{$pair_No} / $nloc_total * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of all loci: missing ind 2";
for $pair_No (sort {$a <=> $b} keys %{$counts{'NN_ind2'}}) {
	$val = $counts{'NN_ind2'}{$pair_No} / $nloc_total * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
print SUMOUT "% of all loci: missing both";
for $pair_No (sort {$a <=> $b} keys %{$counts{'NN_both'}}) {
	$val = $counts{'NN_both'}{$pair_No} / $nloc_total * 100;
	print SUMOUT "\t$val";
}
print SUMOUT "\n";
close SUMOUT;
#}
