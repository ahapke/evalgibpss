#!/usr/bin/perl -w
#comp_simres version 06 Copyright 2015 Andreas Hapke
#This program compares genotypes simulated with
#sim_dat to those inferred by GIbPSs from simulated data
#out of sim_dat or subsamp_e_N_01.pl.
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
use Math::Round qw(nearest);

#Keyword Infilecolumns! marks where column order of infiles matters

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
M => '',#1 analyze merged sequences, 0 analyze unmerged sequences
i => '',#input directory for files about simulated data
o => ''#output directory
);
my %defaults = (
M => '0',#1 analyze merged sequences, 0 analyze unmerged sequences
i => 'sim_dat',#input directory for files about simulated data
o => 'comp_simres'#output directory
);
my $n_arg = 0;#number command line arguments
my $flag = '';#a flag
my $val = '';#a value
#Input settings
my $ori_popall_fname = 'ori_popall.txt';#popall outfile out of sim_dat
my $inf_popall_fname = 'popall.txt';#popall outfile out of poploc
my $ori_popallcol = 0;#column of popallseq in ori_popall.txt, changes according to -M
my $inf_popallcol = 0;#column of popallseq in popall.txt, changes according to -M
my $ori_gt_fname = 'ori_gt.txt';#genotypes outfile out of sim_dat
my $ori_poploc_fname = 'ori_poploc.txt';#out of sim_dat, contains information about indel loci
my $selind_fname = 'export/sel_ind.txt';#individual selection file out of data_selector
my $selgt_fname = 'export/sel_gt.txt';#genotype selection file out of data_selector
my $selloc_fname = 'export/sel_loc.txt';#loci selection file out of data_selector
my $ind_fname = 'individuals.txt';#individuals file out of indloc
my $indpoplocsuff = '_indpoploc.txt';#suffix of indpoploc file out of indpoploc
my $indpopallsuff = '_indpopall.txt';#suffix of indpopall file out of indpoploc
my $inslash = "\/";
my $outslash = "\/";
#Output settings
my $loccount_fname = 'loc_counts.txt';#summary counts over loci
my $indcount_fname = 'ind_gt_counts.txt';#summary counts for individuals
my $percgtinf_fname = 'perc_inf_gt.txt';#percent of genotypes in inferred dataset
#Other variables
my %ori_popall = ();# popallseq = (poplocID,popall_ID)
my %ori_poploc = ();# {poplocID}{popall_ID} = popallseq
my %ori_indel = ();# {poplocID} = 0/1 1: has indel variation
my %ori_dup = ();# {ori_poplocID}=0/1 1: is a duplicate locus
my %inf_popall = ();# popallseq = (poplocID,popall_ID)
my %inf_poploc = ();# {poplocID}{popall_ID} = popallseq
my $ori_popallseq = '';
my $ori_poplocID = 0;
my $ori_popall_ID = 0;
my $inf_popallseq = '';
my $inf_poplocID = 0;
my $inf_popall_ID = 0;
my $for_found = 0;#found forward version
my $rc_found = 0;#found revcomp version
my %ori_poploc_cat = ();# {poplocID} = ex/mis/pres
						#for exclude/missing/present in inf data
my %inf_poploc_cat = ();# {poplocID} = ex/err/corr
						#for exclude/erroneous/correct
my $cat = '';# a category
my %poplocID_oriinf = ();# {ori_poplocID}=inf_poplocID
my %poplocID_infori = ();# {inf_poplocID}=ori_poplocID
my %popall_ID_infori = ();# {inf_popall_ID} = ori_popall_ID / -1 for incorrect
my %ori_popall_thisloc = ();# {ori_popallseq} = ori_popall_ID
my %ori_gt = ();# {ind}{ori_poplocID}{ori_popall_ID} = 1 true genotypes
my %inf_gt = ();# {ori_poplocID}{ori_popall_ID} = 1 inferred genotypes of one ind
				#translated into ori IDs, 
my %sel_ind = ();# {ind}=1
my %sel_gt = ();# {inf_poplocID}{ind}=1
my $ind = '';#an individual ID
my $infilename = '';
my %errloc = ();# {inf_poplocID}=1
my %exloc = ();# {inf_poplocID}=1
my %indcounts = ();#counts from genotype comparison:
	# {ind}{id}=identical
	#      {drop}=allelic dropout
	#      {mm}=mismatch
	#      {ms}=missing
	#      {exp}=excluded poploc
	#      {erp}=erroneous poploc
	#      {misxrp}= excluded and erroneous loci missing in this ind
my %indinfperc = ();#counts in %indcounts as percent of selected loci in %sel_gt
my %mmainfperc = ();# {'min'} {'max'} {'avg'} across inds of percentages in %indinfperc
my %loccounts = ();#counts of loci
    # {ori}: keys ex indelex mis indelmis pres indelpres total totindel
	# {inf}: keys ex err corr total
my $n_selloc = 0;#number of selected loci from inferred data
my $n_ind = 0;#number of selected individuals
my $total = 0;#sum of percentage of all inds
my $id = 0;# a counter
my $tempstring1 = '';
my @temparr1 = ();
my @temparr2 = ();
my %temphash1 = ();
my $i = 0;
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
#M integer 0/1
unless(($user_settings{'M'} eq '0') or ($user_settings{'M'} eq '1')) {
	$user_settings{'M'} = $defaults{'M'};
}
#i, if defined check if it exists, else use default check if exists
if (length $user_settings{'i'} > 0) {
	#remove trailing slash or backslash if any
	if ($user_settings{'i'} =~ m/\//) {
		$user_settings{'i'} =~ s/\/$//;
	}
	elsif ($user_settings{'i'} =~ m/\\/) {
		$user_settings{'i'} =~ s/\\$//;
	}
	if (length $user_settings{'i'} == 0) {
		print "unusable indirname, using default $defaults{'i'}\n";
		$user_settings{'i'} = $defaults{'i'};
	}
	unless (-d "$user_settings{'i'}") {
		print "Directory $user_settings{'i'} does not exist, using default $defaults{'i'}\n";
		$user_settings{'i'} = $defaults{'i'};
	}
} else {
	$user_settings{'i'} = $defaults{'i'};
}
unless(-d $user_settings{'i'}) {
	print "Directory $user_settings{'i'} does not exist, exiting..\n";
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
#}

##############################
#Include pathes into filenames
##############################

#{
if ($user_settings{'i'} =~ m/\\/) {
	$inslash = "\\";
}
if ($user_settings{'o'} =~ m/\\/) {
	$outslash = "\\";
}
$ori_popall_fname = $user_settings{'i'} . $inslash . $ori_popall_fname;
$ori_gt_fname = $user_settings{'i'} . $inslash . $ori_gt_fname;
$ori_poploc_fname = $user_settings{'i'} . $inslash . $ori_poploc_fname;

$loccount_fname = $user_settings{'o'} . $outslash . $loccount_fname;
$indcount_fname = $user_settings{'o'} . $outslash . $indcount_fname;
$percgtinf_fname = $user_settings{'o'} . $outslash . $percgtinf_fname;
#}

#######################################################
#Read in existing preselection or create full selection
#######################################################

#{
#If all three selection files are there
if ((-f $selind_fname) and (-f $selloc_fname) and (-f $selgt_fname)) {
	print "Loading your selection of individuals...\n";
	unless(open(INFILE,$selind_fname)) {#open it
		print "Cannot open $selind_fname. Exiting..\n";
		exit;
	}
	@temparr1 = <INFILE>;#copy file into array
	close INFILE;
	for $ind (@temparr1) {
		chomp $ind;
		$sel_ind{$ind} = 1;#store selected inds in %sel_ind
	}
	#if the selection file was empty (%sel_ind still empty now)
	if (keys %sel_ind == 0) {
		print "Your selection evaluates to an empty dataset.\n",
		"I have nothing to compare Exiting..\n";
		exit;			
	}
	print "Loading your selection of genotypes...\n";
	unless(open(INFILE,$selgt_fname)) {#open it
		print "Cannot open $selgt_fname. Exiting..\n";
		exit;
	}
	while ($tempstring1 = <INFILE>) {#read in line by line
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);#split into array
		$sel_gt{$temparr1[0]}{$temparr1[1]} = 1;#store in %sel_gt
	}
	close INFILE;
	#if the selection file was empty (%sel_gt still empty now)
	if (keys %sel_gt == 0) {
		print "Your selection evaluates to an empty dataset.\n",
		"I have nothing to compare. Exiting..\n";
		exit;			
	}	
} else {#None or not all selection files there
	print "I did not find three selection files in directory export.\n",
	"I will create a complete selection of all inferred data.\n";
	#Create sel_ind from file individuals.txt
	unless(open(INFILE,$ind_fname)) {#open it
		print "Cannot open $ind_fname. Exiting..\n";
		exit;
	}
	@temparr1 = <INFILE>;#copy into array
	close INFILE;
	for $ind (@temparr1) {
		chomp $ind;
		$sel_ind{$ind} = 1;#store in %sel_ind
	}
	#Create sel_gt from *_indpoploc.txt files
	for $ind (keys %sel_ind) {
		#read in file *_indpoploc.txt for this ind (outfile from indpoploc)
		$infilename = $ind . $indpoplocsuff;
		unless(open(INFILE,$infilename)) {#open or die
			print "Cannot open $infilename. Exiting..\n";
			exit;
		}
		$tempstring1 = <INFILE>;#read to skip headerline
		while ($tempstring1 = <INFILE>) {#read data
			chomp $tempstring1;
			@temparr1 = split(/\t/,$tempstring1);
			#Infilecolumns!
			$sel_gt{$temparr1[1]}{$ind} = 1;
		}
		close INFILE;
	}
}
#}

###########################################
#Link poplocs, identify nonmatching poplocs
###########################################

#{
#Analyze merged or unmerged popall sequences
#Infilecolumns!
if ($user_settings{'M'} == 0) {#unmerged
	$ori_popallcol = 3;
	$inf_popallcol = 4;
} else {#merged
	$ori_popallcol = 2;
	$inf_popallcol = 2;
}
#Open ori_popall.txt and popall.txt
unless(open(ORIPOPALL,$ori_popall_fname)) {
	print "Cannot open $ori_popall_fname, exiting..\n";
	exit;
}
unless(open(INFPOPALL,$inf_popall_fname)) {
	print "Cannot open $inf_popall_fname, exiting..\n";
	exit;
}
#Read in ori_popall.txt, populate %ori_poploc {poplocID}{popall_ID}=popallseq
#later: populate %ori_popall: popallseq = (poplocID,popall_ID)
$tempstring1 = <ORIPOPALL>;#skip headerline
while ($tempstring1 = <ORIPOPALL>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	#@{$ori_popall{$temparr1[$ori_popallcol]}} = ($temparr1[0],$temparr1[1]);
	$ori_poploc{$temparr1[0]}{$temparr1[1]} = $temparr1[$ori_popallcol];
}
close ORIPOPALL;
#Read in popall.txt,
#populate %inf_popall: popallseq = (poplocID,popall_ID)
#populate %inf_poploc: {poplocID}{popall_ID} = popallseq
#Exclude all alleles of loci that are not in selection
#Read in all other alleles

$tempstring1 = <INFPOPALL>;#skip headerline
while ($tempstring1 = <INFPOPALL>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	if (defined $sel_gt{$temparr1[0]}) {#if this inf poploc is selected
		@{$inf_popall{$temparr1[$inf_popallcol]}} = ($temparr1[0],$temparr1[1]);
		$inf_poploc{$temparr1[0]}{$temparr1[1]} = $temparr1[$inf_popallcol];
	}
}
close INFPOPALL;
#match ori poplocs to inf poplocs, search forward and reverse complement
#you don't know, which version of the locus is in inferred data
%temphash1 = ();
$for_found = 0;
$rc_found = 0;
for $ori_poplocID (keys %ori_poploc) {
	for $ori_popall_ID (keys %{$ori_poploc{$ori_poplocID}}) {
		$ori_popallseq = $ori_poploc{$ori_poplocID}{$ori_popall_ID};
		#search forward version
		if (defined $inf_popall{$ori_popallseq}) {#found this allele in inferred data (inf)
			$temphash1{$inf_popall{$ori_popallseq}[0]} = 1;#collect inf_poplocID
			$for_found = 1;#found an allele in forward orientation
		}
		#search reverse complement
		$ori_popallseq = reverse $ori_popallseq;
		$ori_popallseq =~ tr/ACGTacgt/TGCAtgca/;
		if (defined $inf_popall{$ori_popallseq}) {#found this allele in inferred data (inf)
			$temphash1{$inf_popall{$ori_popallseq}[0]} = 1;#collect inf_poplocID
			$rc_found = 1;#found an allele in forward orientation
		}			
	}
	if (($for_found + $rc_found) < 2) {#found alleles in one orientation or didn't find alleles
		if (keys %temphash1 == 1) {#one inf locus matches
			$inf_poplocID = (sort keys %temphash1)[0];
			$ori_poploc_cat{$ori_poplocID} = 'pres';#flag ori_poploc as present in inf
			$poplocID_oriinf{$ori_poplocID} = $inf_poplocID;#store matching inf poplocID
		}
		elsif (keys %temphash1 > 1) {#several inf loci match
			$ori_poploc_cat{$ori_poplocID} = 'ex';#flag ori_poploc as exclude
			for $inf_poplocID (keys %temphash1) {
				$inf_poploc_cat{$inf_poplocID} = 'ex';#flag each inf_poploc as exclude
			}
		} else {#no inferred locus matches
			$ori_poploc_cat{$ori_poplocID} = 'mis';#flag ori_poploc as missing in inf
		}
		if ($rc_found == 1) {#found reverse complement
			#change all ori alleles to reverse complement
			for $ori_popall_ID (keys %{$ori_poploc{$ori_poplocID}}) {
				$ori_poploc{$ori_poplocID}{$ori_popall_ID} = reverse $ori_poploc{$ori_poplocID}{$ori_popall_ID};
				$ori_poploc{$ori_poplocID}{$ori_popall_ID} =~ tr/ACGTacgt/TGCAtgca/;
			}
		}
		#populate %ori_popall for this locus popallseq = (poplocID,popall_ID)
		for $ori_popall_ID (keys %{$ori_poploc{$ori_poplocID}}) {
			$ori_popallseq = $ori_poploc{$ori_poplocID}{$ori_popall_ID};
			@{$ori_popall{$ori_popallseq}} = ($ori_poplocID,$ori_popall_ID);
		}
	} else {#found alleles in both orientations
		$ori_poploc_cat{$ori_poplocID} = 'ex';#flag ori_poploc as exclude
		for $inf_poplocID (keys %temphash1) {
			$inf_poploc_cat{$inf_poplocID} = 'ex';#flag each inf_poploc as exclude
		}
		#populate %ori_popall for this locus popallseq = (poplocID,popall_ID)
		#use forward and reverse complement of popallseq
		for $ori_popall_ID (keys %{$ori_poploc{$ori_poplocID}}) {
			$ori_popallseq = $ori_poploc{$ori_poplocID}{$ori_popall_ID};
			@{$ori_popall{$ori_popallseq}} = ($ori_poplocID,$ori_popall_ID);
			$ori_popallseq = reverse $ori_popallseq;
			$ori_popallseq =~ tr/ACGTacgt/TGCAtgca/;
			@{$ori_popall{$ori_popallseq}} = ($ori_poplocID,$ori_popall_ID);
		}
			
	}
	$for_found = 0;
	$rc_found = 0;
	%temphash1 = ();#set back and next ori_poploc
}
#match inf poplocs to ori poplocs
%temphash1 = ();
for $inf_poplocID (keys %inf_poploc) {
	for $inf_popall_ID (keys %{$inf_poploc{$inf_poplocID}}) {
		$inf_popallseq = $inf_poploc{$inf_poplocID}{$inf_popall_ID};
		if (defined $ori_popall{$inf_popallseq}) {#found this allele in ori data
			$temphash1{$ori_popall{$inf_popallseq}[0]} = 1;#collect ori_poplocID
		}
	}
	if (keys %temphash1 == 1) {#one locus matches
		$ori_poplocID = (sort keys %temphash1)[0];
		#unless this inf poploc is already flagged as exclude
		unless((defined $inf_poploc_cat{$inf_poplocID}) and ($inf_poploc_cat{$inf_poplocID} eq 'ex')) {
			$inf_poploc_cat{$inf_poplocID} = 'corr';#flag as correct locus
			$poplocID_infori{$inf_poplocID} = $ori_poplocID;#store matching ori_poplocID
		}
	}
	elsif (keys %temphash1 > 1) {#several ori loci match
		$inf_poploc_cat{$inf_poplocID} = 'ex';#flag inf poploc as exclude
		for $ori_poplocID (keys %temphash1) {
			$ori_poploc_cat{$ori_poplocID} = 'ex';#flag each ori_poploc as exclude
		}
	} else {#no ori locus matches
		$inf_poploc_cat{$inf_poplocID} = 'err';#flag inf poploc as erroneous
	}
	%temphash1 = ();#set back and next inf_poploc
}

#Check %poplocID_oriinf again and remove loci now flagged as exclude
for $ori_poplocID (keys %poplocID_oriinf) {
	if ($ori_poploc_cat{$ori_poplocID} eq 'ex') {
		delete $poplocID_oriinf{$ori_poplocID};
	}
}
#}

###############################################################
#Read in ori_poploc.txt
#populate %ori_indel {ori_poplocID}=0/1, 1: has indel variation
#populate %ori_dup {ori_poplocID}=0/1 1: is a duplicate locus
###############################################################

#{
unless(open(ORIPL, $ori_poploc_fname)) {
	print "Cannot open $ori_poploc_fname, exiting..\n";
	exit;
}
$tempstring1 = <ORIPL>;
while ($tempstring1 = <ORIPL>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	#Infilecolumns!
	$ori_indel{$temparr1[0]} = $temparr1[1];
	$ori_dup{$temparr1[0]} = $temparr1[2];
}
close ORIPL;
#}

#######################
#Determine locus counts
#######################

#{
$loccounts{'ori'}{'totindel'} = 0;
$loccounts{'ori'}{'totdup'} = 0;
$loccounts{'ori'}{'ex'} = 0;
$loccounts{'ori'}{'indelex'} = 0;
$loccounts{'ori'}{'dupex'} = 0;
$loccounts{'ori'}{'mis'} = 0;
$loccounts{'ori'}{'indelmis'} = 0;
$loccounts{'ori'}{'dupmis'} = 0;
$loccounts{'ori'}{'pres'} = 0;
$loccounts{'ori'}{'indelpres'} = 0;
$loccounts{'ori'}{'duppres'} = 0;
$loccounts{'inf'}{'ex'} = 0;
$loccounts{'inf'}{'err'} = 0;
$loccounts{'inf'}{'corr'} = 0;
for $ori_poplocID (keys %ori_poploc_cat) {
	$cat = $ori_poploc_cat{$ori_poplocID};
	if ($ori_indel{$ori_poplocID} == 1) {
		$cat = 'indel' . $cat;
		++$loccounts{'ori'}{'totindel'};
	}
	if ($ori_dup{$ori_poplocID} == 1) {
		$cat = 'dup' . $cat;
		++$loccounts{'ori'}{'totdup'};
	}
	++$loccounts{'ori'}{$cat};
}
for $inf_poplocID (keys %inf_poploc_cat) {
	$cat = $inf_poploc_cat{$inf_poplocID};
	++$loccounts{'inf'}{$cat};
}
$loccounts{'ori'}{'total'} = keys %ori_poploc_cat;
$loccounts{'inf'}{'total'} = keys %inf_poploc_cat;
#}

#######################################
#Link popalls
#Match inf_popallseqs to ori_popallseqs
#locus-wise, for matching loci only
#######################################

#{
for $inf_poplocID (keys %poplocID_infori) {
	$ori_poplocID = $poplocID_infori{$inf_poplocID};
	#collect all ori_popallseqs of this locus
	for $ori_popall_ID (keys %{$ori_poploc{$ori_poplocID}}) {
		$ori_popallseq = $ori_poploc{$ori_poplocID}{$ori_popall_ID};
		$ori_popall_thisloc{$ori_popallseq} = $ori_popall_ID;
	}
	#look up inf_popallseqs in %ori_popall_thisloc
	for $inf_popall_ID (keys %{$inf_poploc{$inf_poplocID}}) {
		$inf_popallseq = $inf_poploc{$inf_poplocID}{$inf_popall_ID};
		if (defined $ori_popall_thisloc{$inf_popallseq}) {
			$ori_popall_ID = $ori_popall_thisloc{$inf_popallseq};
			$popall_ID_infori{$inf_popall_ID} = $ori_popall_ID;
		} else {
			$popall_ID_infori{$inf_popall_ID} = -1;
		}
	}	
	%ori_popall_thisloc = ();#set back and next inf_poploc
}
#}

#########################################
#Read in ori_gt.txt
#Exclude excluded and missing ori_poplocs
#########################################

#{
unless(open(ORIGT,$ori_gt_fname)) {
	print "Cannot open $ori_gt_fname, exiting..\n";
	exit;
}
$tempstring1 = <ORIGT>;
while ($tempstring1 = <ORIGT>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	if ($ori_poploc_cat{$temparr1[1]} eq 'pres') {
		#Infilecolumns!
		$ori_gt{$temparr1[0]}{$temparr1[1]}{$temparr1[2]} = 1;
		$ori_gt{$temparr1[0]}{$temparr1[1]}{$temparr1[3]} = 1;
	}
}
#}

##############################
#Compare genotypes of each ind
##############################

#{
for $ind (sort keys %sel_ind) {

	#read in file *_indpopall, populate %inf_gt {ori_poplocID}{ori_popall_ID} = 1
	$infilename = $ind . $indpopallsuff;
	unless(open(INDPOPALL,$infilename)) {
		print "Cannot open $infilename, exiting..\n";
		exit;
	}
	$tempstring1 = <INDPOPALL>;
	while ($tempstring1 = <INDPOPALL>) {
		chomp $tempstring1;
		@temparr1 = split(/\t/,$tempstring1);
		#Infilecolumns!
		#If locus is selected for this ind
		if ((defined $sel_gt{$temparr1[2]}) and (defined $sel_gt{$temparr1[2]}{$ind})) {
			$cat = $inf_poploc_cat{$temparr1[2]};#cat: ex/err/corr
			if ($cat eq 'ex') {#if locus is excluded
				$exloc{$temparr1[2]} = 1;#collect inf_poplocID to count later
			}
			elsif ($cat eq 'err') {#if locus is erroneous
				$errloc{$temparr1[2]} = 1;#collect inf_poplocID to count later
			} else {#locus is correct
				#Translate inf_poplocID and inf_popall_ID into ori IDs
				#ori_popall_ID can be -1 when inf_popall is incorrect
				$ori_poplocID = $poplocID_infori{$temparr1[2]};
				$ori_popall_ID = $popall_ID_infori{$temparr1[3]};
				$inf_gt{$ori_poplocID}{$ori_popall_ID} = 1;#store genotype
			}
		}
	}
	close INDPOPALL;
	
	#Count excluded and erroneous loci in this ind
	$indcounts{$ind}{'exp'} = keys %exloc;
	$indcounts{$ind}{'erp'} = keys %errloc;
	%errloc = ();#set back
	%exloc = ();#set back
	
	#Compare genotypes
	$indcounts{$ind}{'id'} = 0;
	$indcounts{$ind}{'drop'} = 0;
	$indcounts{$ind}{'mm'} = 0;
	$indcounts{$ind}{'ms'} = 0;
	for $ori_poplocID (keys %{$ori_gt{$ind}}) {
		if (defined $inf_gt{$ori_poplocID}) {#locus found in inf data
			if (defined $inf_gt{$ori_poplocID}{-1}) {#if any allele in inferred data incorrect
				++$indcounts{$ind}{'mm'};#count mismatch
			}
			#number of ori alleles 2, inf alleles: 1 =>check for dropout or mismatch
			elsif((keys %{$ori_gt{$ind}{$ori_poplocID}} == 2) and
			(keys %{$inf_gt{$ori_poplocID}} == 1)) {
				$ori_popall_ID = (sort keys %{$inf_gt{$ori_poplocID}})[0];
				if (defined $ori_gt{$ind}{$ori_poplocID}{$ori_popall_ID}) {
					++$indcounts{$ind}{'drop'};#count allelic dropout
				} else {
					++$indcounts{$ind}{'mm'};#count mismatch
				}
			}
			#number of alleles different not 2 alleles in ori 1 allele in inf
			elsif (keys %{$ori_gt{$ind}{$ori_poplocID}} != keys %{$inf_gt{$ori_poplocID}}) {
				++$indcounts{$ind}{'mm'};#count mismatch
			} else {#compare all alleles
				$id = 1;
				for $ori_popall_ID (keys %{$ori_gt{$ind}{$ori_poplocID}}) {
					unless(defined $inf_gt{$ori_poplocID}{$ori_popall_ID}) {
						$id = 0;
					}
				}
				if ($id == 0) {
					++$indcounts{$ind}{'mm'};#count mismatch
				} else {
					++$indcounts{$ind}{'id'};#count identical genotype
				}
			}
		} else {#locus not found in inf data
			++$indcounts{$ind}{'ms'};#count missing genotype
		}
	}
	#Determine number of excluded and erroneous loci missing in this ind
	$indcounts{$ind}{'misxrp'} = $loccounts{'inf'}{'ex'} + $loccounts{'inf'}{'err'} -
	$indcounts{$ind}{'exp'} - $indcounts{$ind}{'erp'};
	%inf_gt = ();#set back and next ind
}
#}

####################################################
#Determine percentages of genotypes in inferred data
#Determine min, max, average percentages across inds
####################################################

#{
$n_selloc = keys %sel_gt;#number of selected loci: 100%
#individual percentages, not rounded
for $ind (keys %indcounts) {
	for $cat (keys %{$indcounts{$ind}}) {
		$indinfperc{$ind}{$cat} = $indcounts{$ind}{$cat} / $n_selloc * 100;
	}
}
$n_ind = keys %sel_ind;
#collect cats in %indinfperc
@temparr1 = ();
$ind = (sort keys %indinfperc)[0];
for $cat (keys %{$indinfperc{$ind}}) {
	push @temparr1, $cat;
}
#min max mean of individual percentages, rounded
$total = 0;
@temparr2 = ();
for $cat (@temparr1) {
	for $ind (keys %indinfperc) {
		$total += $indinfperc{$ind}{$cat};
		push @temparr2, $indinfperc{$ind}{$cat};
	}
	@temparr2 = sort {$a <=> $b} @temparr2;
	$mmainfperc{'min'}{$cat} = nearest(0.001,$temparr2[0]);
	$mmainfperc{'max'}{$cat} = nearest(0.001,pop @temparr2);
	$mmainfperc{'avg'}{$cat} = nearest(0.001,($total / $n_ind));
	$total = 0;
	@temparr2 = ();
}
@temparr1 = ();
#Round individual percentages
for $ind (keys %indcounts) {
	for $cat (keys %{$indcounts{$ind}}) {
		$indinfperc{$ind}{$cat} = nearest(0.001,$indinfperc{$ind}{$cat});
	}
}
#}

#################
#Produce outfiles
#################

#{
#Locus counts
unless(open(OUTLC, ">$loccount_fname")) {
	print "Cannot open $loccount_fname, exiting..\n";
	exit;
}
print OUTLC "original loci total\t$loccounts{'ori'}{'total'}\n";
print OUTLC "original loci with indel\t$loccounts{'ori'}{'totindel'}\n";
print OUTLC "original loci with duplicate\t$loccounts{'ori'}{'totdup'}\n";
print OUTLC "found original loci without indel\t$loccounts{'ori'}{'pres'}\n";
print OUTLC "found original loci with indel\t$loccounts{'ori'}{'indelpres'}\n";
print OUTLC "found original loci with duplicate\t$loccounts{'ori'}{'duppres'}\n";
print OUTLC "excluded original loci without indel\t$loccounts{'ori'}{'ex'}\n";
print OUTLC "excluded original loci with indel\t$loccounts{'ori'}{'indelex'}\n";
print OUTLC "excluded original loci with duplicate\t$loccounts{'ori'}{'dupex'}\n";
print OUTLC "missing original loci without indel\t$loccounts{'ori'}{'mis'}\n";
print OUTLC "missing original loci with indel\t$loccounts{'ori'}{'indelmis'}\n";
print OUTLC "missing original loci with duplicate\t$loccounts{'ori'}{'dupmis'}\n";
print OUTLC "inferred loci\t$loccounts{'inf'}{'total'}\n";
print OUTLC "correctly inferred loci\t$loccounts{'inf'}{'corr'}\n";
print OUTLC "excluded inferred loci\t$loccounts{'inf'}{'ex'}\n";
print OUTLC "erroneous inferred loci\t$loccounts{'inf'}{'err'}\n";
close OUTLC;

#Individual genotype counts
unless(open(OUTIC, ">$indcount_fname")) {
	print "Cannot open $indcount_fname, exiting..\n";
	exit;
}
print OUTIC "ind\tidentical\tdropout\tmismatch\tmissing\texcluded_loc\terr_loc\tmissing_excl_and_err\n";
for $ind (sort keys %indcounts) {
	print OUTIC "$ind",
	"\t$indcounts{$ind}{'id'}",
	"\t$indcounts{$ind}{'drop'}",
	"\t$indcounts{$ind}{'mm'}",
	"\t$indcounts{$ind}{'ms'}",
	"\t$indcounts{$ind}{'exp'}",
	"\t$indcounts{$ind}{'erp'}",
	"\t$indcounts{$ind}{'misxrp'}\n";
}
close OUTIC;

#Percent of genotypes in inferred dataset
unless(open(OUTPC, ">$percgtinf_fname")) {
	print "Cannot open $percgtinf_fname, exiting..\n";
	exit;
}
print OUTPC "ind\tidentical\tdropout\tmismatch\tmissing\texcluded_loc\terr_loc\tmissing_excl_and_err\n";
@temparr1 = ('min','max','avg');
for $ind (@temparr1) {
	print OUTPC "$ind",
	"\t$mmainfperc{$ind}{'id'}",
	"\t$mmainfperc{$ind}{'drop'}",
	"\t$mmainfperc{$ind}{'mm'}",
	"\t$mmainfperc{$ind}{'ms'}",
	"\t$mmainfperc{$ind}{'exp'}",
	"\t$mmainfperc{$ind}{'erp'}",
	"\t$mmainfperc{$ind}{'misxrp'}\n";
}
for $ind (sort keys %indinfperc) {
	print OUTPC "$ind",
	"\t$indinfperc{$ind}{'id'}",
	"\t$indinfperc{$ind}{'drop'}",
	"\t$indinfperc{$ind}{'mm'}",
	"\t$indinfperc{$ind}{'ms'}",
	"\t$indinfperc{$ind}{'exp'}",
	"\t$indinfperc{$ind}{'erp'}",
	"\t$indinfperc{$ind}{'misxrp'}\n";
}
close OUTPC;
#}
#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";

exit;
