#!/usr/bin/perl -w
#sim_fastq version 02 Copyright 2016 Andreas Hapke
#This program simulates sequence data that can be analyzed with GIbPSs and pyRAD.
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
use List::Util qw(shuffle);
use IO::File; #requires Perl 5.004 or higher
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);


#uses subroutines
#simloc
#simall
#simfreq
#simread
#getfastqfn
#checkphredlen
#generate_bc
#permute_bc
#swap
#bc_ok
#dist

#runtime parameters
my ($start_time) = time();
my ($end_time) = 0;
my ($run_s) = 0;

#######################
#Declare and initialize
#######################

#{
#Variables for user settings
my %user_settings = (
nI => '',#number of individuals
nL => '',#number of loci
nR => '',#number of reads per locus and ind
Lsl => '',#locus sequence length
Ld => '',#minimum distance between loci
Rl => '',#read length (default: int(Lsl/2)
nM => '',#number of mutation steps for creation of alleles
e => '',#polymerase error rate
rc => '',#Boolean: print reads as revcomp with prob. 0.5
fqi => '',#input directory with fastq files
o => '',#output directory
nsf => '',#maximum number of sequences per fastq outfile
rs => ''#random number seed
);
my %defaults = (
nI => '10',#number of individuals
nL => '100',#number of loci
nR => '40',#number of reads per locus and ind
Lsl => '120',#locus sequence length
Ld => '0',#minimum distance between loci
nM => '0',#number of mutation steps for creation of alleles
e => '0.001',#polymerase error rate
rc => '0',#Boolean: print reads as revcomp with prob. 0.5
o => 'sim_fastq',#output directory
nsf => '1000000',#maximum number of sequences per fastq outfile
rs => ''#random number seed
);
my $n_arg = 0;#number command line arguments
my $flag = '';#a flag
my $val = '';#a value
#Output settings
my $repfilename = 'sim_fastq_rep.txt';#Name of report file
my $popall_fname = 'ori_popall.txt';#Name of popall outfile
my $poploc_fname = 'ori_poploc.txt';#Name of poploc outfile
my $gt_fname = 'ori_gt.txt';#Name of genotypes outfile
my $popallfas_fname = 'popall.fas';#Name of fasta file with all alleles
my $inds_bc_fname = 'inds_barcodes.txt';#Name of file with ind IDs and barcodes for pyRAD
my $bc_fname = 'barcodes.txt';#Name of file with barcodes for GIbPSs (fdm)
my $seq_err_fname = 'seq_err_rate.txt';#Name of file with sequencing error rates
my $slash = "\/";
#Other variables
my $seed = '';
my @fastqfinfn = ();#names of fastq f-infiles;
my @fastqrinfn = ();#names of fastq r-infiles;
my %popall = ();# {poplocID}{popall_ID} = popallseq
my %popallfreq = ();# {poplocID}{popall_ID} = popallfreq
					#popallfreq: count in population of nI inds
my %poploc_indel = ();# {poplocID} = 0/1 1 for has indel variation
my %poploc_dup = ();# {poplocID} = 0/1 1 for locus has a duplicate
my %gt = ();# {ind}{poplocID}=(popall_ID,popall_ID)
my %bcs = ();# {barcode}=1 barcodes
my %ind_bc = ();# {ind} = barcode
my %fseqerr = ();# {pos} = error count, pos starts with 0
my %rseqerr = ();# {pos} = error count, pos starts with 0
my $noutseq = 0;#Total number of sequence pairs to produce
my $pos = 0;#read position
my $outpos = 0;#pos+1
my $ferate = 0;#error rate in f
my $rerate = 0;#error rate in r
my $bclen = 6;#barcode length
my $bc = '';#barcode
my $ind = '';#ID of an individual
my $poplocID = 0;#ID of a locus
my $popall_ID = 0;#ID of an allele
my $popallseq = '';#sequence of an allele
my $popallseqnm = '';#unmerged sequence of an allele
my $frestr = 'GATCC';#forward restriction site, must have same length as $rrestr
my $rrestr = 'GGATC';#reverse restriction site, must have same length as $frestr
my $outreadlen = 0;#final read length with barcode, restriction site, and one additional char (pyRAD may eat one char)
my $f_read = '';#forward read
my $r_read = '';#reverse read
my $r_start = 0;#start position of r_read (in f-strand orientation)
my $head = '';#fastq header
my $fqname = '';#name of a fastq file
my $i = 0;
my $j = 0;
my @temparr1 = ();
#}

########################
#Take over user-settings
########################

#{
if (@ARGV) {#if user provided some
	$n_arg = @ARGV;#determine number of arguments
	for ($i = 0; $i < ($n_arg - 1); $i += 2) {
		$flag = $ARGV[$i];
		$flag =~ s/^-//;
		$val = $ARGV[$i+1];
		#if flag is defined, take over,
		if (defined $user_settings{$flag}) {
			$user_settings{$flag} = $val;
		}
	}
}
#Check user settings for plausibility
#nI integer > 0
unless(($user_settings{'nI'} =~ /^\d+$/) and ($user_settings{'nI'} > 0)) {
	$user_settings{'nI'} = $defaults{'nI'};
}
#nL integer > 0
unless(($user_settings{'nL'} =~ /^\d+$/) and ($user_settings{'nL'} > 0)) {
	$user_settings{'nL'} = $defaults{'nL'};
}
#nR integer > 0
unless(($user_settings{'nR'} =~ /^\d+$/) and ($user_settings{'nR'} > 0)) {
	$user_settings{'nR'} = $defaults{'nR'};
}
#Lsl integer >= 20
unless(($user_settings{'Lsl'} =~ /^\d+$/) and ($user_settings{'Lsl'} >= 20)) {
	$user_settings{'Lsl'} = $defaults{'Lsl'};
}
#Ld integer >=0, <= Lsl/3
unless(($user_settings{'Ld'} =~ /^\d+$/) and ($user_settings{'Ld'} >= 0) and
($user_settings{'Ld'} <= ($user_settings{'Lsl'} / 3))) {
	$user_settings{'Ld'} = $defaults{'Ld'};
}
#Rl integer >= 10, <= Lsl, default int(Lsl/2)
unless(($user_settings{'Rl'} =~ /^\d+$/) and ($user_settings{'Rl'} >= 10) and
($user_settings{'Rl'} <= $user_settings{'Lsl'})) {
	$user_settings{'Rl'} = int($user_settings{'Lsl'} / 2);
}
#nM integer >= 0, <= Lsl
unless(($user_settings{'nM'} =~ /^\d+$/) and ($user_settings{'nM'} >= 0) and
($user_settings{'nM'} <= $user_settings{'Lsl'})) {
	$user_settings{'nM'} = $defaults{'nM'};
}
#e decimal >= 0, <= 1
unless((($user_settings{'e'} =~ /^\d$/) or ($user_settings{'e'} =~ /^\d\.\d+$/))
and ($user_settings{'e'} >= 0) and ($user_settings{'e'} <= 1)) {
	$user_settings{'e'} = $defaults{'e'};
}
#rc must be 0 or 1
unless(($user_settings{'rc'} eq '0') or ($user_settings{'rc'} eq '1')) {
	$user_settings{'rc'} = $defaults{'rc'};
}
#fqi must be defined and must exist
unless(length $user_settings{'fqi'} > 0){
	print "Directory with fastq infiles must be defined (flag -fqi), exiting..\n";
	exit;
}
if($user_settings{'fqi'} =~ /\\/){
	$user_settings{'fqi'} =~ s/\\$//;#remove trailing backslash
}
elsif($user_settings{'fqi'} =~ /\//){
	$user_settings{'fqi'} =~ s/\/$//;#remove trailing slash
}
unless(-d $user_settings{'fqi'}){
	print "Cannot find directory $user_settings{'fqi'}, exiting..\n";
	exit;
}
#o if defined, try to create or exit; if not defined, try to create default output dir or exit
if (length $user_settings{'o'} > 0) {
	#remove trailing slash or backslash if any
	if ($user_settings{'o'} =~ m/\//) {
		$user_settings{'o'} =~ s/\/$//;
	}
	elsif ($user_settings{'o'} =~ m/\\/) {
		$user_settings{'o'} =~ s/\\$//;
	}
	if (length $user_settings{'o'} == 0) {
		print "unusable outdirname, using default $defaults{'o'}\n";
		$user_settings{'o'} = $defaults{'o'};
	}
} else {
	$user_settings{'o'} = $defaults{'o'};
}
if (-d "$user_settings{'o'}")  {
	print "Directory $user_settings{'o'} already exists.\n",
	"Please rename or delete. Exiting..\n";
	exit;
} else {
	unless (mkdir "$user_settings{'o'}") {
		print "Cannot create directory $user_settings{'o'}, exiting..\n";
		exit;
	}
}
#nsf must be integer > 0
if(length $user_settings{'nsf'} > 0){
	unless(($user_settings{'nsf'} =~ /^\d+$/) and ($user_settings{'nsf'} > 0)){
		$user_settings{'nsf'} = $defaults{'nsf'};
		print "Max number of sequences per fastq outfile (-nsf) must be integer > 0,\n",
		"using default: $user_settings{'nsf'}.\n";
	}
} else {
	$user_settings{'nsf'} = $defaults{'nsf'};
	print "Using default max number of sequences per fastq outfile: $user_settings{'nsf'}.\n";
}
#Check that not more than 1000 fastq outfile pairs will be produced
$i = $user_settings{'nI'} * $user_settings{'nL'} * $user_settings{'nR'};#number of sequences
if(($i % $user_settings{'nsf'}) > 0){
	$j = int($i / $user_settings{'nsf'}) + 1;
} else{
	$j = $i / $user_settings{'nsf'};
}
if($j > 1000){
	print "I cannot produce more than 1000 pairs of fastq outfiles, exiting..\n";
	exit;
}
#rs must be positive integer
if (length $user_settings{'rs'} > 0) {
	unless (($user_settings{'rs'} =~ /^\d+$/) and ($user_settings{'rs'} > 0)) {
		print "Setting -rs: must be a positive integer, exiting..\n";
		exit;
	}
}
#}

#######################
#Get random number seed
#######################

#{
$seed = $user_settings{'rs'};
if (length $seed > 0) {
	srand($seed);
} else {
	$seed = srand();
}
#}

#######################################
#Include outdirectory into outfilenames
#######################################

#{
if ($user_settings{'o'} =~ m/\\/) {#if outdirname contains backslash
	$slash = "\\";
}
$repfilename = $user_settings{'o'} . $slash . $repfilename;
$popall_fname = $user_settings{'o'} . $slash . $popall_fname;
$poploc_fname = $user_settings{'o'} . $slash . $poploc_fname;
$gt_fname = $user_settings{'o'} . $slash . $gt_fname;
$popallfas_fname = $user_settings{'o'} . $slash . $popallfas_fname;
$inds_bc_fname = $user_settings{'o'} . $slash . $inds_bc_fname;
$bc_fname = $user_settings{'o'} . $slash . $bc_fname;
$seq_err_fname = $user_settings{'o'} . $slash . $seq_err_fname;
#}

################################################################################
#Get fastq infilenames, determine $outreadlen,
#Check if Phred lines in fastq-infiles have sufficient length (must be all equal
################################################################################

#{
getfastqfn($user_settings{'fqi'},\@fastqfinfn,\@fastqrinfn);
$outreadlen = $user_settings{'Rl'} + (length $frestr) + $bclen + 10;
checkphredlen($outreadlen,\@fastqfinfn,\@fastqrinfn);
#}

#######################################################
#Open report file and print settings to file and screen
#Print fastq infilenames to report file
#######################################################

#{
unless(open(OUTREP, ">", $repfilename)) {
	print "Cannot open outfile $repfilename, exiting ...\n\n";
	exit;
}
binmode OUTREP;
print OUTREP "SIM_FASTQ used settings:\12";
print OUTREP "nI $user_settings{'nI'}\12";
print OUTREP "nL $user_settings{'nL'}\12";
print OUTREP "nR $user_settings{'nR'}\12";
print OUTREP "Lsl $user_settings{'Lsl'}\12";
print OUTREP "Ld $user_settings{'Ld'}\12";
print OUTREP "Rl $user_settings{'Rl'}\12";
print OUTREP "nM $user_settings{'nM'}\12";
print OUTREP "e $user_settings{'e'}\12";
print OUTREP "fqi $user_settings{'fqi'}\12";
print OUTREP "o $user_settings{'o'}\12";
print OUTREP "nsf $user_settings{'nsf'}\12";
print OUTREP "rc $user_settings{'rc'}\12";
print OUTREP "rs $seed\12";

print OUTREP "\12Forward fastq infiles:\12";
for $fqname (@fastqfinfn){
	print OUTREP "$fqname\12";
}
print OUTREP "\12Reverse fastq infiles:\12";
for $fqname (@fastqrinfn){
	print OUTREP "$fqname\12";
}

print  "Used settings:\n";
print  "nI $user_settings{'nI'}\n";
print  "nL $user_settings{'nL'}\n";
print  "nR $user_settings{'nR'}\n";
print  "Lsl $user_settings{'Lsl'}\n";
print  "Ld $user_settings{'Ld'}\n";
print  "Rl $user_settings{'Rl'}\n";
print  "nM $user_settings{'nM'}\n";
print  "e $user_settings{'e'}\n";
print  "fqi $user_settings{'fqi'}\n";
print  "o $user_settings{'o'}\n";
print  "nsf $user_settings{'nsf'}\n";
print  "rc $user_settings{'rc'}\n";
print  "rs $seed\n";
#}

##################
#Generate barcodes
##################

#{
generate_bc($user_settings{'nI'},\%bcs,$frestr);
$i = keys %bcs;
unless($i >= $user_settings{'nI'}){
	print "I found only $i barcodes",
	" - cannot simulate $user_settings{'nI'} individuals, exiting..\n";
	exit;
}
#}

##############
#Open outfiles
##############

#{
unless(open(OUTPOPALL, ">", $popall_fname)) {
	print "Cannot open outfile $popall_fname, exiting ..\n";
	exit;
}
binmode OUTPOPALL;
unless(open(OUTPOPLOC, ">", $poploc_fname)) {
	print "Cannot open outfile $poploc_fname, exiting..\n";
	exit;
}
binmode OUTPOPLOC;
unless(open(OUTGT, ">", $gt_fname)) {
	print "Cannot open outfile $gt_fname, exiting..\n";
	exit;
}
binmode OUTGT;
unless(open(OUTFAS, ">", $popallfas_fname)) {
	print "Cannot open outfile $popallfas_fname, exiting ...\n\n";
	exit;
}
binmode OUTFAS;
unless(open(OUTINDSBC, ">", $inds_bc_fname)) {
	print "Cannot open outfile $inds_bc_fname, exiting ...\n\n";
	exit;
}
binmode OUTINDSBC;
unless(open(OUTBC, ">", $bc_fname)) {
	print "Cannot open outfile $bc_fname, exiting ...\n\n";
	exit;
}
binmode OUTBC;
unless(open(OUTSR, ">", $seq_err_fname)) {
	print "Cannot open outfile $seq_err_fname, exiting ...\n\n";
	exit;
}
binmode OUTSR;
#}

#############################
#Simulate data,print fq files
#############################

#{
#Create first allele of each locus, populate %popall {poplocID}{popall_ID} = popallseq
print "Simulating loci..\n";
simloc(\%popall,$user_settings{'nL'},$user_settings{'Lsl'},
$user_settings{'Ld'},$user_settings{'nM'});
#Create nM additional alleles for each locus, add to %popall
print "Simulating alleles..\n";
simall(\%popall,$user_settings{'Lsl'},$user_settings{'nM'});
#Set "locus has duplicate" to 0 for all poplocs
for $poplocID (keys %popall) {
	$poploc_dup{$poplocID} = 0;
}
#Create random allele frequencies for alleles (including 0)
#populate %popallfreq {poplocID}{popall_ID} = popallfreq
simfreq(\%popall,\%popallfreq,$user_settings{'nI'});
#Remove alleles with frequency 0
for $poplocID (keys %popallfreq) {
	for $popall_ID (keys %{$popallfreq{$poplocID}}) {
		if ($popallfreq{$poplocID}{$popall_ID} == 0) {
			delete $popall{$poplocID}{$popall_ID};
			delete $popallfreq{$poplocID}{$popall_ID};
		}
	}
}
#Set "has indel" to 0 for all poplocs
for $poplocID (keys %popall) {
	$poploc_indel{$poplocID} = 0;
}
#Sample genotypes
print "Simulating genotypes..\n";
simgt(\%popallfreq,$user_settings{'nI'},\%gt);

#Assign barcodes to individuals
@temparr1 = (sort keys %bcs);
for $ind (sort keys %gt){
	$bc = shift @temparr1;
	$ind_bc{$ind} = $bc;
}

#Create reads and print to fastq files
print "Simulating reads..\n";
simread(\%popall,\%gt,$user_settings{'nR'},$user_settings{'Rl'},$outreadlen,
$user_settings{'e'},$user_settings{'rc'},$user_settings{'o'},$slash,
$user_settings{'nsf'},$frestr,$rrestr,\%ind_bc,\@fastqfinfn,\@fastqrinfn,
\%fseqerr,\%rseqerr);
#}

########################
#Print to other outfiles
########################

#{
#Print data to popall outfile
print OUTPOPALL "poplocID\tpopall_ID\tpopallseq_merged\tpopallseq_notmerged\tfreq\12";
for $poplocID (sort {$a <=> $b} keys %popall) {
	for $popall_ID (sort {$a <=> $b} keys %{$popall{$poplocID}}) {
		$popallseq = $popall{$poplocID}{$popall_ID};
		$f_read = substr($popallseq,0,$user_settings{'Rl'});
		$r_start = (length $popallseq) - $user_settings{'Rl'};
		$r_read = substr($popallseq,$r_start,$user_settings{'Rl'});
		$popallseqnm = $f_read . $r_read;
		print OUTPOPALL "$poplocID\t$popall_ID\t$popallseq\t$popallseqnm\t$popallfreq{$poplocID}{$popall_ID}\12";
	}
}
close OUTPOPALL;
#print data to poploc outfile
print OUTPOPLOC "poplocID\tindel\tdup\12";
for $poplocID (sort {$a <=> $b} keys %poploc_indel) {
	print OUTPOPLOC "$poplocID\t$poploc_indel{$poplocID}\t$poploc_dup{$poplocID}\12";
}
close OUTPOPLOC;
#print data to gt outfile
print OUTGT "ind\tpoplocID\tpopall_1\tpopall_2\12";
for $ind (sort keys %gt) {
	for $poplocID (sort {$a <=> $b} keys %{$gt{$ind}}) {
		print OUTGT "$ind\t$poplocID\t$gt{$ind}{$poplocID}[0]\t$gt{$ind}{$poplocID}[1]\12";
	}
}
close OUTGT;
#print data to popallfas outfile
for $poplocID (sort {$a <=> $b} keys %popall) {
	for $popall_ID (sort {$a <=> $b} keys %{$popall{$poplocID}}) {
		$popallseq = $popall{$poplocID}{$popall_ID};
		$head = '>' . $poplocID . '_' . $popall_ID;
		print OUTFAS "$head\12";
		print OUTFAS "$popallseq\12";
	}
}
#print data to inds barcodes outfile and barcodes outfile
for $ind (sort keys %ind_bc){
	print OUTINDSBC "$ind\t$ind_bc{$ind}\12";
	print OUTBC "$ind_bc{$ind}\12";
}
close OUTINDSBC;
close OUTBC;
#print data to seqerr outfile
print OUTSR "pos\tf\tr\12";
$noutseq = $user_settings{'nI'} * $user_settings{'nL'} * $user_settings{'nR'};
for $pos(sort {$a <=> $b} keys %fseqerr){
	$outpos = $pos + 1;
	$ferate = $fseqerr{$pos} / $noutseq;
	$rerate = $rseqerr{$pos} / $noutseq;
	print OUTSR "$outpos\t$ferate\t$rerate\12";
}
close OUTSR;
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
print OUTREP "Analysis completed.\12Run took $run_s seconds.\12\12";
close OUTREP;

exit;

############
#Subroutines
############

#definition of subroutine simloc
#Creates random loci (one allele per locus)
#Loci have a configurable minimum distance to each other
#all have same length and equal base compositions (as close as possible)
#Expects:
#ref to %popall {poplocID}{popall_ID} = seq		will hold IDs and sequences
#nL  number of loci
#Lsl sequence length
#Ld  minimum distance between two loci
#nM  number of additional alleles that will be created for each locus later
#
#poplocIDs start with 1 and increment by 1
#popall_IDs start with 1 and increment by nM+1
sub simloc {
	#Declare and initialize, _r means a reference
	my ($popall_r,$nL,$Lsl,$Ld,$nM) = @_;
	my $poplocID = 1;
	my $popall_ID = 1;
	my $nA = 0;#number of A in sequence
	my $nC = 0;#number of C in sequence
	my $nG = 0;#number of G in sequence
	my $nT = 0;#number of T in sequence
	my $startseq = '';#starting sequence
	my @seq = ();#current sequence as array
	my $currseq = '';#current sequence of a locus
	my %poploc = ();# {poplocID}=seq
	my $dist = 0;#number of differences between two seqs
	my $i = 0;

	#create startseq
	$nA = int($Lsl / 4);
	$nC = int($Lsl / 4);
	$nG = int($Lsl / 4);
	$nT = $Lsl - $nA - $nC -$nG;
	$startseq = 'A' x $nA;
	$startseq .= 'C' x $nC;
	$startseq .= 'G' x $nG;
	$startseq .= 'T' x $nT;
	@seq = split('', $startseq);
	@seq = shuffle @seq;

	#create loci
	$dist = $Ld;
	while ($poplocID <= $nL) {
		@seq = shuffle @seq;#randomize characters in seq
		$currseq = join('',@seq);
		#test if $currseq is sufficiently distant from all known ones
		for $poplocID (keys %poploc) {
			$dist = $Lsl - (($currseq ^ $poploc{$poplocID}) =~ tr/\0//);
			last if ($dist < $Ld);
		}
		if ($dist >= $Ld) {#if sequence was sufficiently distant
			$poploc{$poplocID} = $currseq;
			++$poplocID;
		}
	}
	#populate %popall owned by main
	for $poplocID (sort {$a <=> $b} keys %poploc) {
		$$popall_r{$poplocID}{$popall_ID} = $poploc{$poplocID};
		$popall_ID += ($nM + 1);
	}	
}

#definition of subroutine simall
#Takes a set of loci, each with one allele (seq: ACGT, all capital)
#For each locus: creates nM additional alleles
#by stepwise mutation: one mutation per mutation step, at random position
#ensures that all alleles are different
#nM: number of mutation steps
#Expects:
#ref to %popall  {poplocID}{popall_ID}=seq contains one allele per locus
#Lsl sequence length (equal for all loci)
#nM  number of mutation steps
#adds the new alleles to %popall
#poplocIDs are ascending from 1
#popall_IDs in input start with 1 and increase by nM+1
#new popall IDs increment by 1 from that of first allele of each locus
sub simall {
	#Declare and initialize, _r means a reference
	my ($popall_r,$Lsl,$nM) = @_;
	my %subst = ();# {character}=array of substition characters, e.g. {A}=(C,G,T)
	my $ori_char = '';#original character at a position
	my $sub_char = '';#substition character
	my $poplocID = 0;
	my $popall_ID = 0;
	my $popallseq = '';
	my $mutpos = 0;#position of a mutation
	my $i = 0;
	my $j = 0;
	my %alleles = ();# {seq} = 1 alleles of one locus
	my $unique = 0;# 1 if a new allele is different from all others of the locus
	
	#initialize %subst
	@{$subst{'A'}} = ('C','G','T');
	@{$subst{'C'}} = ('A','G','T');
	@{$subst{'G'}} = ('A','C','T');
	@{$subst{'T'}} = ('A','C','G');
	
	for $poplocID (sort {$a <=> $b} keys %{$popall_r}) {#loop through loci
		$popall_ID = (sort {$a <=> $b} keys %{$$popall_r{$poplocID}})[0];
		$popallseq = $$popall_r{$poplocID}{$popall_ID};#get first allele
		$alleles{$popallseq} = 1;#store in %alleles
		for ($i = 0; $i < $nM; ++$i) {#perform $nM steps
			$unique = 0;
			while ($unique == 0) {
				$mutpos = int(rand $Lsl);#select a random position
				$ori_char = substr($popallseq,$mutpos,1);
				$j = int(rand 3);#select a random number 0..2 to choose a substition character
				$sub_char = $subst{$ori_char}[$j];#get substitution character
				substr($popallseq,$mutpos,1) = $sub_char;#substitute
				#test if new allele is different from already existing ones
				unless(defined $alleles{$popallseq}) {
					$unique = 1;
					++$popall_ID;
					$$popall_r{$poplocID}{$popall_ID} = $popallseq;#store new allele
					$alleles{$popallseq} = 1;
				}
			}
		}
		%alleles = ();
	}
}

#definition of subroutine simfreq
#Creates random allele frequencies for each allele
#Random allele frequencies include 0
#Frequencies are counts in population of nI individuals
#If a locus has at least 2 alleles, ensures that at least
#2 alleles get frequency > 0
#Expects
#ref to %popall {poplocID}{popall_ID} = popallseq
#ref to %popallfreq {poplocID}{popall_ID} = popallfreq
#$nI  number of individuals
#populates %popallfreq
sub simfreq {
	#Declare and initialize, _r means a reference
	my ($popall_r, $popallfreq_r, $nI) = @_;
	my $poplocID = 0;
	my $popall_ID = 0;
	my @all_IDs = ();#all popall_IDs of a locus
	my $first_all_ID = 0;#ID of first allele of a locus after order randomization
	my $first_all_freq = 0;#frequency of first allele
	my $last_all_ID = 0;#ID of last allele of a locus after order randomization
	my $last_all_freq = 0;#frequency of last allele
	my $freq = 0;#frequency
	my $remcount = 0;#total remaining count for alleles that have not yet a frequency
	
	for $poplocID (sort {$a <=> $b} keys %{$popall_r}) {
		$remcount = 2 * $nI - 1;
		$last_all_freq = 1;
		@all_IDs = (sort {$a <=> $b} keys %{$$popall_r{$poplocID}});#get popall_IDs of this loc
		if (@all_IDs > 1) {#if there is more than one allele, make sure that at least 2 alleles have frequency > 0 
			@all_IDs = shuffle @all_IDs;#randomize order of alleles of this loc
			$first_all_ID = shift @all_IDs;#separate first allele
			$first_all_freq = 1;#set frequency of first allele to 1
			$remcount -= 1;#reduce remaining count accordingly
			$freq = int(rand ($remcount + 1));#determine a random frequency
			$first_all_freq += $freq;#add to frequency of first allele
			$remcount -= $freq;#reduce remaining count accordingly
			$$popallfreq_r{$poplocID}{$first_all_ID} = $first_all_freq;
		}
		$last_all_ID = pop @all_IDs;#separate last allele
		for $popall_ID (@all_IDs) {
			$freq = int(rand ($remcount + 1));
			$$popallfreq_r{$poplocID}{$popall_ID} = $freq;
			$remcount -= $freq;
		}
		$last_all_freq += $remcount;
		$$popallfreq_r{$poplocID}{$last_all_ID} = $last_all_freq;
	}	
}

#definition of subroutine simgt
#Creates genotypes
#Expects:
#ref to %popallfreq {poplocID}{popall_ID}=freq (count of allele in nI individuals)
#$nI number of individuals
#ref to %gt {ind}{poplocID}=(popall_ID,popall_ID)
#
#Sum of freq for one poploc in %popallfreq must be equal to 2*$nI
#Populates %gt
sub simgt {
	#Declare and initialize, _r means a reference
	my ($popallfreq_r,$nI,$gt_r) = @_;
	my $poplocID = 0;
	my $popall_ID = 0;
	my $popall_ID1 = 0;
	my $popall_ID2 = 0;
	my @gt = ();
	my $freq = 0;
	my @popalls = ();
	my %inds = ();# {ind} = 1
	my $ind = '';
	my $i = 0;
	
	#Create individual IDs
	for ($i = 1; $i <= $nI; ++$i) {
		$ind = sprintf('%04d',$i);
		$ind = 'ind' . $ind;
		$inds{$ind} = 1;
	}
	#Create genotypes
	for $poplocID (sort {$a <=> $b} keys %{$popallfreq_r}) {
		for $popall_ID (sort {$a <=> $b} keys %{$$popallfreq_r{$poplocID}}) {
			$freq = $$popallfreq_r{$poplocID}{$popall_ID};
			#build array of popall_IDs, each one $freq times in array
			for ($i = 0; $i < $freq; ++$i) {
				push @popalls, $popall_ID;
			}
		}
		@popalls = shuffle @popalls;#randomize order of popall_IDs
		#draw individual genotypes for this loc
		for $ind (sort keys %inds) {
			$popall_ID1 = shift @popalls;
			$popall_ID2 = shift @popalls;
			@gt = sort {$a <=> $b} ($popall_ID1,$popall_ID2);
			@{$$gt_r{$ind}{$poplocID}} = @gt;
		}		
	}
}

#Definition of subroutine simread
#Simulates paired reads and prints them to paired fastq files
#For each ind and locus:
#Takes a genotype with two allele sequences
#Creates nR allele sequences, with random proportions of the two alleles
#optionally takes each sequence as revcomp with prob. 0.5
#appends restriction sites to both ends
#adds barcode at beginning
#simulates polymerase errors with probability e
#simulates f- and r-reads with length outreadlen
#gets a pair of real Phred score lines from fastq infiles
#simulates sequencing errors based on Phred scores with offset 33
#counts sequencing errors (from Phred scores) and populates 
#%fseqerr and %rseqerr owned by main
#Expects
#ref to %popall {poplocID}{popall_ID}=popallseq
#ref to %gt {ind}{poplocID}=(popall_ID,popall_ID) popall_IDs sorted ascending
#$nR	number of reads per locus and ind
#$Rl:	final read length
#$Rlbcr	read length with barcode and restriction sites
#$e		polymerase error rate
#$rc    1: print reads as reverse complement with prob. 0.5, 0: don't
#$outdirname name of outfile directory for fastq files
#$slash contains slash or backslash for building of fastq outfilenames
#$nsf	maximum number of sequences per fastq outfile
#$frestr forward restriction site
#$rrestr reverse restriction site
#ref to %ind_bc {ind}=barcode
#ref to @fastqfinfn fastq f-infilenames
#ref to @fastqrinfn fastq r-infilenames
#ref to %fseqerr {pos} = error count, pos starts with 0
#ref to %rseqerr {pos} = error count, pos starts with 0
sub simread {
	#Declare and initialize, _r means a reference
	my ($popall_r,$gt_r,$nR,$Rl,$Rlbcr,$e,$rc,$outdirname,$slash,
	$nsf,$frestr,$rrestr,$ind_bc_r,$fastqfinfn_r,$fastqrinfn_r,
	$fseqerr_r,$rseqerr_r) = @_;
	my $ind = '';
	my $poplocID = 0;
	my $popallID = 0;
	my $popallseq = '';
	my $f_read = '';
	my $r_start = 0;
	my $r_read = '';
	my $readNo = 0;
	my $f_head = '';
	my $r_head = '';
	my $nind = 0;#number of individuals
	my $nloc = 0;#number of loci
	my $nseqmax = 0;#total number of sequences (f-only) to produce
	my $nseqnow = 0;#total number of sequences (f-only) produced so far
	my %pord_pe = ();# {Phred score + offset}=error prob
	my $poffset = 33;#Phred score offset for Sanger
	my $phred = 0;#Phred score
	my $pord = 0;#order number (Phred score plus offset)
	my $pe = 0;#Phred error probability
	my $Nord = 0;#order number that corresponds to Phred 2
	my $foutfn = '';#fastq-f-outfilename
	my $routfn = '';#fastq-r-outfilename
	my $outfno = 1;#fastq outfile number
	my $outseqno = 0;#number of sequences in current fastq outfiles
	my $finfn = '';#fastq-f-infilename
	my $finfh = '';#fastq-f-infilehandle
	my $rinfn = '';#fastq-r-outfilename
	my $rinfh = '';#fastq-r-outfilehandle
	my $infno = 0;#fastq infile number
	my @gt = ();#one genotype (two allele IDs)
	my @alleles = ();#both allele IDs of genotype, each nR times
	my $pos = 0;#position in sequence
	my %error = ();# {ori_char}=array of possible err_chars
	my $ori_char = '';#original character at a position
	my $err_char = '';#erroneous character
	my $finline = '';
	my $rinline = '';
	my $fpline = '';#Phred score line
	my $rpline = '';#Phred score line
	my $i = 0;
	my $j = 0;
	
	#Determine $nseqmax
	$nind = keys %{$gt_r};
	$nloc = keys %{$popall_r};
	$nseqmax = $nind * $nloc * $nR;
	#Initialize %error
	@{$error{'A'}} = ('C','G','T');
	@{$error{'C'}} = ('A','G','T');
	@{$error{'G'}} = ('A','C','T');
	@{$error{'T'}} = ('A','C','G');
	#Initialize sequencing error counts
	for($pos = 0; $pos < $Rlbcr; ++$pos){
		$$fseqerr_r{$pos} = 0;
		$$rseqerr_r{$pos} = 0;
	}
	#populate %pord_pe and $Nord
	for ($phred = 3; $phred <= 50; ++$phred){
		$pe = 10 ** ($phred / -10);
		$pord_pe{($phred + $poffset)} = $pe;
	}
	$Nord = 2 + $poffset;	
	#open first pair of fastq outfiles
	$foutfn = sprintf('%04d',$outfno);
	++$outfno;
	$foutfn .= '.fastq';
	$foutfn = 'sfq_R1_' . $foutfn;
	$routfn = $foutfn;
	$routfn =~ s/_R1_/_R2_/;
	$foutfn = $outdirname . $slash . $foutfn;
	$routfn = $outdirname . $slash . $routfn;
	unless(open(FQFOUT,">",$foutfn)){
		print "Cannot open outfile $foutfn, exiting..\n";
		exit;
	}
	unless(open(FQROUT,">",$routfn)){
		print "Cannot open outfile $routfn, exiting..\n";
		exit;
	}
	binmode FQFOUT;#File will be written in binary mode, enables Linux newlines
	binmode FQROUT;#File will be written in binary mode, enables Linux newlines
	#open first pair of fastq infiles and read first line from each
	$finfn = $$fastqfinfn_r[$infno];
	$rinfn = $$fastqrinfn_r[$infno];
	++$infno;
	if($infno == @{$fastqfinfn_r}){#reached last fastq infile
		$infno = 0;#restart with first file next time
	}
	if($finfn =~ /.gz/){
		$finfh = new IO::Uncompress::Gunzip $finfn;
	} else {
		$finfh = IO::File->new("< $finfn");
	}
	if($rinfn =~ /.gz/){
		$rinfh = new IO::Uncompress::Gunzip $rinfn;
	} else {
		$rinfh = IO::File->new("< $rinfn");
	}
	unless(defined $finfh){
		print "Cannot open infile $finfn, exiting..\n";
		exit;
	}
	unless(defined $rinfh){
		print "Cannot open infile $rinfn, exiting..\n";
		exit;
	}
	$finline = <$finfh>;
	$rinline = <$rinfh>;
	#create reads
	for $ind (sort keys %{$gt_r}){
		print "Individual $ind\n";
		for $poplocID (sort {$a <=> $b} keys %{$$gt_r{$ind}}) {
			@gt = @{$$gt_r{$ind}{$poplocID}};
			#sample each popall_ID of genotype nR times into @alleles
			for ($i = 0; $i < $nR; ++$i) {
				push @alleles, $gt[0];
				push @alleles, $gt[1];
			}
			@alleles = shuffle @alleles;#randomize their order
			for ($i = 0; $i < $nR; ++$i) {
				$popallseq = $$popall_r{$poplocID}{$alleles[$i]};
				#if reads shall be saved as revcomp with probability 0.5,do
				if (($rc == 1) and ((rand) < 0.5)) {
					$popallseq = reverse $popallseq;
					$popallseq =~ tr/ACGTacgt/TGCAtgca/;
				}
				#Append barcode and restriction site
				$popallseq = $$ind_bc_r{$ind} . $frestr . $popallseq . $rrestr;
				#If polymerase error rate > 0: Simulate polymerase errors
				if ($e > 0) {
					for ($pos = 0; $pos < length $popallseq; ++$pos) {
						if ((rand) < $e) {#select error position with probability $e
							$ori_char = substr($popallseq,$pos,1);
							$j = int(rand 3);#select a random number 0..2 to choose an error character
							$err_char = $error{$ori_char}[$j];#get error character
							substr($popallseq,$pos,1) = $err_char;#substitute
						}
					}
				}
				#Construct reads
				$f_read = substr($popallseq,0,$Rlbcr);
				$r_start = (length $popallseq) - $Rlbcr;
				$r_read = substr($popallseq,$r_start,$Rlbcr);
				$r_read = reverse $r_read;
				$r_read =~ tr/ACGTacgt/TGCAtgca/;
				#Get Phred lines and truncate to $Rlbcr (we have already read the header line)
				$finline = <$finfh>;
				$finline = <$finfh>;
				$finline = <$finfh>;
				$finline =~ s/(\n|\r\n?)//;
				$fpline = substr($finline,0,$Rlbcr);
				$rinline = <$rinfh>;
				$rinline = <$rinfh>;
				$rinline = <$rinfh>;
				$rinline =~ s/(\n|\r\n?)//;
				$rpline = substr($rinline,0,$Rlbcr);
				#Simulate sequencing errors in f_read
				for($pos = 0; $pos < length $fpline; ++$pos){
					$pord = ord(substr($fpline,$pos,1));#includes offset
					if($pord == $Nord){#Phred score is 2
						substr($f_read,$pos,1) = 'N';
					} else {
						$pe = $pord_pe{$pord};
						if((rand) < $pe){
							$ori_char = substr($f_read,$pos,1);
							$j = int(rand 3);#select a random number 0..2 to choose an error character
							$err_char = $error{$ori_char}[$j];#get error character
							substr($f_read,$pos,1) = $err_char;#substitute
							++$$fseqerr_r{$pos};#count error
						}
					}
				}
				#Simulate sequencing errors in r_read
				for($pos = 0; $pos < length $rpline; ++$pos){
					$pord = ord(substr($rpline,$pos,1));#includes offset
					if($pord == $Nord){#Phred score is 2
						substr($r_read,$pos,1) = 'N';
					} else {
						$pe = $pord_pe{$pord};
						if((rand) < $pe){
							$ori_char = substr($r_read,$pos,1);
							$j = int(rand 3);#select a random number 0..2 to choose an error character
							$err_char = $error{$ori_char}[$j];#get error character
							substr($r_read,$pos,1) = $err_char;#substitute
							++$$rseqerr_r{$pos};#count error
						}
					}
				}
				#create f- and r-header
				$readNo = $i + 1;
				$f_head = '@HWI-ST558:200:C3BVRACXX:4:' .
				$poplocID . ':' . $alleles[$i] . ':' . $readNo;
				$r_head = $f_head;
				$f_head .= ' 1:N:0:ATCACG';
				$r_head .= ' 2:N:0:ATCACG';
				#print reads to fastq-outfiles
				print FQFOUT "$f_head\12$f_read\12+\12$fpline\12";
				print FQROUT "$r_head\12$r_read\12+\12$rpline\12";
				++$outseqno;
				++$nseqnow;
				#open new fastq outfiles if necessary
				if($outseqno == $nsf and $nseqnow < $nseqmax){
					$outseqno = 0;
					close FQFOUT;
					close FQROUT;
					$foutfn = sprintf('%04d',$outfno);
					++$outfno;
					$foutfn .= '.fastq';
					$foutfn = 'sfq_R1_' . $foutfn;
					$routfn = $foutfn;
					$routfn =~ s/_R1_/_R2_/;
					$foutfn = $outdirname . $slash . $foutfn;
					$routfn = $outdirname . $slash . $routfn;
					unless(open(FQFOUT,">",$foutfn)){
						print "Cannot open outfile $foutfn, exiting..\n";
						exit;
					}
					unless(open(FQROUT,">",$routfn)){
						print "Cannot open outfile $routfn, exiting..\n";
						exit;
					}
					binmode FQFOUT;#File will be written in binary mode, enables Linux newlines
					binmode FQROUT;#File will be written in binary mode, enables Linux newlines
				}
				#read first lines from current or next pair of fastq infiles
				unless(defined <$finfh> and defined <$rinfh>){
					close $finfh;
					close $rinfh;
					#open next pair of fastq infiles and read first line from each
					$finfn = $$fastqfinfn_r[$infno];
					$rinfn = $$fastqrinfn_r[$infno];
					++$infno;
					if($infno == @{$fastqfinfn_r}){#reached last fastq infile
						$infno = 0;#restart with first file next time
					}
					if($finfn =~ /.gz/){
						$finfh = new IO::Uncompress::Gunzip $finfn;
					} else {
						$finfh = IO::File->new("< $finfn");
					}
					if($rinfn =~ /.gz/){
						$rinfh = new IO::Uncompress::Gunzip $rinfn;
					} else {
						$rinfh = IO::File->new("< $rinfn");
					}
					unless(defined $finfh){
						print "Cannot open infile $finfn, exiting..\n";
						exit;
					}
					unless(defined $rinfh){
						print "Cannot open infile $rinfn, exiting..\n";
						exit;
					}
					$finline = <$finfh>;
					$rinline = <$rinfh>;
				}
			}
			@alleles = ();#set back
		}
	}
	close FQFOUT;
	close FQROUT;
	close $finfh;
	close $rinfh;
}

#Definition of subroutine getfastqfn
#Expects:
#$indirname name of directory with fastq infiles
#ref to @fastqfinfn names of fastq f-infiles
#ref to @fastqrinfn names of fastq r-infiles
#Searches files with extension fastq or fastq.gz in indir
#populates the two arrays
#Checks if there are corresponding pairs of fastq infiles
sub getfastqfn {
	#Declare and initialize, _r means a reference
	my($indirname,$fastqfinfn_r,$fastqrinfn_r) = @_;
	my $filename = '';
	my $nf = '';
	my $nr = '';
	my $slash = '/';
	my @temparr1 = ();
	my $i = 0;

	if($indirname =~ /\\/){
	$slash = '\\';
	}
	$indirname .= $slash;
	@temparr1 = (glob "$indirname*.fastq");
	if(@temparr1 == 0){
		@temparr1 = (glob "$indirname*.fastq.gz");
	}
	if(@temparr1 == 0){
		print "No files with name ending '.fastq' or '.fastq.gz' in directory $indirname, exiting..\n";
		exit;
	}
	for $filename (@temparr1){
		if($filename =~ /_R1_/){
			push @{$fastqfinfn_r}, $filename;
		}
		elsif($filename =~ /_R2_/){
			push @{$fastqrinfn_r}, $filename;
		}
		else {
			print "At least one *.fastq or *.fastq.gz filename in directory ",
			"$indirname does not contain '_R1_' or '_R2_', exiting..\n";
			exit;
		}
	}
	$nf = @{$fastqfinfn_r};
	$nr = @{$fastqrinfn_r};
	if($nf != $nr){
		print "Unequal numbers of f- and r- fastq infiles, exiting..\n";
		exit;
	}
	@{$fastqfinfn_r} = sort @{$fastqfinfn_r};
	@{$fastqrinfn_r} = sort @{$fastqrinfn_r};
	for ($i = 0; $i < @{$fastqfinfn_r}; ++$i){
		$filename = $$fastqfinfn_r[$i];
		$filename =~ s/_R1_/_R2_/;
		unless($filename eq $$fastqrinfn_r[$i]){
			print "Fastq infilenames do not correspond, exiting..\n";
			exit;
		}
	}
}

#Definition of subroutine checkphredlen
#Expects
#$outreadlen final read length
#ref to @fastqfinfn fastqf-infilenames
#ref to @fastqrinfn fastqr-infilenames
#Checks if first sequence in each fastq-infile
#has length >= outreadlen
sub checkphredlen {
	#Declare and initialize, _r means a reference
	my($outreadlen,$fastqfinfn_r,$fastqrinfn_r) = @_;
	my @infiles = ();
	my $infn = '';#infilename
	my $infh = '';#infilehandle
	my $inline = '';
	
	@infiles = @{$fastqfinfn_r};
	for $infn (@{$fastqrinfn_r}){
		push @infiles, $infn;
	}
	for $infn(@infiles){
		if ($infn =~ /.gz$/){
			$infh = new IO::Uncompress::Gunzip $infn;
		} else {
			$infh = IO::File->new("< $infn");
		}
		unless(defined $infh){
			print "Cannot open $infn, exiting..\n";
			exit;
		}
		$inline = <$infh>;
		$inline = <$infh>;
		$inline = <$infh>;
		$inline = <$infh>;
		$inline =~ s/(\n|\r\n?)//;
		unless((length $inline) >= $outreadlen){
			print "First Phred score line in file $infn too short, exiting..\n";
			exit;
		}
		close $infh;
	}
}

#Definition of sub generate_bc
#generates barcodes of length 6
#with minimum distance of at least 2 to one another
#barcodes have no identical substring of length 4 at equal position
#barcodes have distance of at least 2 to restriction site in any substring
#restriction site may not be longer than 6
#populates %bcs owned by calling sub
#expects
#$nbcmax: maximum number of barcodes to generate
#ref to %bcs: {bc} = 1
#$restr: restriction site
sub generate_bc{
	#Declare and initialize, _r means a reference
	my($nbcmax,$bcs_r,$restr) = @_;
	my @charstrings = (
	'ACGTAC',
	'ACGTGT',
	'ACGTAG',
	'ACGTCT');
	my $charstring = '';
	my @chars = ();
	my $nchars = 0;
	
	if((length $restr) > 6){
		print "Restriction site longer than 6 characters, exiting..\n";
		exit;
	}
	
	for $charstring (@charstrings){
		@chars = split(//,$charstring);#@chars will be modified by sub permute_bc
		$nchars = @chars;
		permute_bc($nchars,\@chars,$bcs_r,$nbcmax,$restr);
		if(keys %{$bcs_r} >= $nbcmax) {
			last;
		}
	}
}

#Definition of sub permute_bc
#Creates barcodes by permutation of an array of characters
#uses Heap's algorithm for permutation
#uses sub bc_ok to check if a barcode is acceptable
#populates a hash of barcodes
#expects
#$n number of characters to permute (first n characters in array)
#ref to @chars: characters to permute
#ref to %bcs: {bc}=1 hash to populate
#$nbcmax: maximum number of barcodes to generate
#$restr: restriction site
sub permute_bc{
	my($n,$chars_r,$bcs_r,$nbcmax,$restr) = @_;
	my $m = $n - 1;
	my $i = 0;
	my $bc = '';
	if($n ==1 ){
		$bc = join("",@{$chars_r});
		if(bc_ok($bc,$bcs_r,$restr)){
		$$bcs_r{$bc} = 1;
		}
	}
	else {
		for($i = 0; $i < $m; ++$i){
			permute_bc($m,$chars_r,$bcs_r,$nbcmax,$restr);
			if(keys %{$bcs_r} >= $nbcmax){
				return;
			}
			if($n%2 == 0){
				swap(\$$chars_r[$i],\$$chars_r[$m]);
			}
			else {
				swap(\$$chars_r[0],\$$chars_r[$m]);
			}
		}
		permute_bc($m,$chars_r,$bcs_r,$nbcmax,$restr);
		if(keys %{$bcs_r} >= $nbcmax){
			return;
		}
	}
}

#Definition of subroutine swap
#swaps values of two scalars that are passed by reference
sub swap{
	my($a_r,$b_r) = @_;
	my $c = 0;
	$c = $$a_r;
	$$a_r = $$b_r;
	$$b_r = $c;
}

#Definition of subroutine bc_ok
#Checks if a candidate barcode is acceptable
#returns 1 if yes and 0 if not
#expects
#$candbc: the candidate barcode
#ref to %bcs: {bc} = 1 hash of existing barcodes (all of equal length as candidate)
#$restr: restriction site
#candidate is acceptable if it
#is new
#has a distance of at least 2 to every existing barcode
#has no identical substring of length 4 with any existing barcode at same position
#has no substring with length of restriction site and distance < 2 to restriction site
sub bc_ok{
	my($candbc,$bcs_r,$restr) = @_;
	my $bc = '';
	my $i = 0;
	if(defined $$bcs_r{$candbc}){
		return 0;
	}
	for $bc (keys %{$bcs_r}){
		if (dist($candbc,$bc) < 2){
			return 0;
		}
		for ($i = 0; $i <= ((length $candbc) - 4); ++$i){
			if(substr($candbc,$i,4) eq substr($bc,$i,4)){
				return 0;
			}
		}
	}
	for($i = 0; $i <= ((length $candbc) - (length $restr)); ++$i){
		if(dist(substr($candbc,$i,length $restr),$restr) < 2){
			return 0;
		}
	}
	return 1;
}

#definition of subroutine dist
#expects 2 strings (must be of equal length for correct result!)
#returns number of differing characters between them
sub dist {
	#declare and initialize
	my ($seq1,$seq2) = @_;
	my $dist = 0;#number of different characters
	
	$dist = (length $seq1) - (($seq1 ^ $seq2) =~ tr/\0//);
	return $dist;	
}













