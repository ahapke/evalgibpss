#!/usr/bin/perl -w
#subsamp_e_N version 01 Copyright 2015 Andreas Hapke
#This program subsamples reads from sequence data simulated with sim_dat.
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
f => '',#name of file with ind TAB name of fastq infile in each line
nRi => '',#number of reads per locus in infiles
nRo => '',#number of reads per locus in outfiles
e => '',#sequencing error rate
N => '',#N rate
o => '',#output directory
rs => ''#random number seed
);
my %defaults = (
f => '',#name of file with ind TAB name of fastq infile in each line
nRi => '',#number of reads per locus in infiles
nRo => '',#number of reads per locus in outfiles
e => '0.001',#sequencing error rate
N => '0',#N rate
o => 'subsamp_e_N',#output directory
rs => ''#random number seed
);
my $n_arg = 0;#number command line arguments
my $flag = '';#a flag
my $val = '';#a value
#Variables for output
my $slash = "\/";
my $repfilename = 'subsamp_e_N_rep.txt';#Name of report file
my $indsinfiles_fname = 'inds_infiles.txt';#Name of inds infiles outfile

#Other variables
my $infilename = '';
my $outfilename = '';
my $seed = 0;
my %inds_infiles = ();# {indID}=fastq infilename
my %inds_outfiles = ();# {indID}=fastq outfilename
my $ind = '';
my $seqline1 = '';
my $seqline2 = '';
my $seqline3 = '';
my $seqline4 = '';
my %error = ();#{ori_char}=array of possible err_chars
@{$error{'A'}} = ('C','G','T');
@{$error{'C'}} = ('A','G','T');
@{$error{'G'}} = ('A','C','T');
@{$error{'T'}} = ('A','C','G');
my $pos = 0;
my $ori_char = '';
my $err_char = '';


my @inlines = ();
my $tempstring1 = '';
my @temparr1 = ();
my $i = 0;
my $j = 0;
my $k = 0;
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
#-f can you open the file?
if (length $user_settings{'f'} > 0) {
	$infilename = $user_settings{'f'};
	unless(open(INFILE, "<", $infilename)) {
		print "I cannot open file $infilename, exiting..\n";
		exit;
	}
	close INFILE;
} else {
	print "You must provide the name of a file about the fastq-infiles with flag -f:\n",
	"Format: one line per infile: indID TAB infilename/path\n",
	"Exiting..\n";
	exit;
}
#-nRi must be positive integer
unless(($user_settings{'nRi'} =~ /^\d+$/) and ($user_settings{'nRi'} > 0)) {
	print "-nRi must be a positive integer, exiting..\n";
	exit;
}
#-nRo must be a positive integer <= nRi
unless(($user_settings{'nRo'} =~ /^\d+$/) and ($user_settings{'nRi'} <= $user_settings{'nRi'})) {
	$user_settings{'nRo'} = $user_settings{'nRi'};
}
#e decimal >= 0, <= 1
unless((($user_settings{'e'} =~ /^\d$/) or ($user_settings{'e'} =~ /^\d\.\d+$/))
and ($user_settings{'e'} >= 0) and ($user_settings{'e'} <= 1)) {
	$user_settings{'e'} = $defaults{'e'};
}
#N decimal >= 0, <= 1
unless((($user_settings{'N'} =~ /^\d$/) or ($user_settings{'N'} =~ /^\d\.\d+$/))
and ($user_settings{'N'} >= 0) and ($user_settings{'N'} <= 1)) {
	$user_settings{'N'} = $defaults{'N'};
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
$indsinfiles_fname = $user_settings{'o'} . $slash . $indsinfiles_fname;
#}

#######################################################
#Open report file and print settings to file and screen
#######################################################

#{
unless(open(OUTREP, ">$repfilename")) {
	print "Cannot open file $repfilename, exiting ...\n\n";
	exit;
}
print OUTREP "subsamp_e_N used settings:\n";
print OUTREP "f $user_settings{'f'}\n";
print OUTREP "nRi $user_settings{'nRi'}\n";
print OUTREP "nRo $user_settings{'nRo'}\n";
print OUTREP "e $user_settings{'e'}\n";
print OUTREP "o $user_settings{'o'}\n";
print OUTREP "rs $seed\n";

print  "subsamp_e_N used settings:\n";
print  "f $user_settings{'f'}\n";
print  "nRi $user_settings{'nRi'}\n";
print  "nRo $user_settings{'nRo'}\n";
print  "e $user_settings{'e'}\n";
print  "o $user_settings{'o'}\n";
print  "rs $seed\n";

#}

###########################
#Read in file about infiles
#Populate %inds_infiles
#Populate %inds_outfiles
###########################

#{
$infilename = $user_settings{'f'};
unless(open(INFILE, "<", $infilename)) {
	print "Cannot open $infilename, exiting..\n";
	exit;
}
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	unless(@temparr1 == 2) {
		print "Unusable line in file $infilename: $tempstring1\n",
		"Exiting..\n";
		print OUTREP "Unusable line in file $infilename:\n$tempstring1\nexited\n";
		exit;
	}
	$ind = $temparr1[0];
	$infilename = $temparr1[1];
	$inds_infiles{$ind} = $infilename;
	$outfilename = $user_settings{'o'} . $slash . $ind . '.fq';
	$inds_outfiles{$ind} = $outfilename;
}
close INFILE;
#}

###################################
#Produce outfile $indsinfiles_fname
###################################

#{
unless(open(OUTFILE, ">", $indsinfiles_fname)) {
	print "Cannot open $indsinfiles_fname, exiting..\n";
	print OUTREP "Cannot open $indsinfiles_fname, exiting..\n";
	exit;
}
for $ind (sort keys %inds_outfiles) {
	print OUTFILE "$ind\t$inds_outfiles{$ind}\n";
}
close OUTFILE;
#}

########################################
#Read in fastq infiles, produce outfiles
########################################

#{
for $ind (sort keys %inds_infiles) {
	$infilename = $inds_infiles{$ind};
	$outfilename = $inds_outfiles{$ind};
	unless(open(FQIN, "<", $infilename)) {
		print "Cannot open $infilename, exiting..\n";
		print OUTREP "Cannot open $infilename, exiting..\n";
		exit;
	}
	unless(open(FQOUT, ">", $outfilename)) {
		print "Cannot open $outfilename, exiting..\n";
		print OUTREP "Cannot open $outfilename, exiting..\n";
		exit;
	}
	while ($tempstring1 = <FQIN>) {
		push @inlines, $tempstring1;
		for($i = 0; $i < (($user_settings{'nRi'} * 4) - 1); ++$i) {#read in all reads of locus
			$tempstring1 = <FQIN>;
			push @inlines, $tempstring1;
		}
		for ($j = 0; $j < $user_settings{'nRo'}; ++$j) {#treat subsample of reads and print to output
			$seqline1 = shift @inlines;
			$seqline2 = shift @inlines;
			$seqline3 = shift @inlines;
			$seqline4 = shift @inlines;
			chomp $seqline2;
			if ($user_settings{'e'} > 0) {#If error rate > 0: Create sequencing error
				for ($pos = 0; $pos < length $seqline2; ++$pos) {
					if ((rand) < $user_settings{'e'}) {#select error position with probability $e
						$ori_char = substr($seqline2,$pos,1);
						$k = int(rand 3);#select a random number 0..2 to choose an error character
						$err_char = $error{$ori_char}[$k];#get error character
						substr($seqline2,$pos,1) = $err_char;#substitute
					}
				}
			}
			#If N rate > 0: Create missing characters
			if ($user_settings{'N'} > 0) {
				for ($pos = 0; $pos < length $seqline2; ++$pos) {
					if ((rand) < $user_settings{'N'}) {#select N position with probability $N
						substr($seqline2,$pos,1) = 'N';#substitute
					}
				}
			}
			#print sequence to outfile
			print FQOUT "$seqline1$seqline2\n$seqline3$seqline4";
		}
		@inlines=();
	}
	close FQIN;
	close FQOUT;
}
#}


#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
print OUTREP "Analysis completed.\nRun took $run_s seconds.\n\n";
close OUTREP;

exit;
