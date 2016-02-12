#!/usr/bin/perl -w
#remtrim version 02 Copyright 2015 Andreas Hapke
#This program reads *.edit files produced by pyRAD version 3.0.1.
#It filters the sequences for trimmed reads and prints them
#to a new directory called editsnotrim.
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


#Declare and initialize
my $inds_fname = 'inds.barcodes';
my @inds = ();#holds individual IDs
my $ind = '';#ind ID
my $editstem = 'edits/';#first part of path to edits files
my $editex = '.edit';#extension of edits files
my $line1 = '';
my $line2 = '';
my $newdirname = 'editsnotrim';
my $editfname = '';
my $neweditfname = '';
my $repfilename = 'trimrep.txt';
my %countstrim = ();# {ind} = count of discarded trimmed reads
my %countsgood = ();# {ind} = count of good retained reads
my $tempstring1 = '';
my @temparr1 = ();



unless(open(RUNTIMEOUT, ">>", "runtime.txt")) {
	print "Cannot open outfile, exiting..\n";
	exit;
}
print RUNTIMEOUT "remtrim_01: remove trimmed reads from edits files\n";


#read in barcodes file
unless(open(INFILE, "<", $inds_fname)) {
	print "Cannot open $inds_fname, exiting..\n";
	exit;
}
while ($tempstring1 = <INFILE>) {
	chomp $tempstring1;
	@temparr1 = split(/\t/,$tempstring1);
	push @inds, $temparr1[0];
}
close INFILE;

#create new directory
if (-d "$newdirname") {
	print "Directory $newdirname already exists.\n",
	"Please remove or rename, exiting..\n";
	exit;
}
unless(mkdir "$newdirname") {
	print "Cannot create directory $newdirname, exiting..\n";
	exit;
}

#initialize counts
for $ind (@inds) {
	$countstrim{$ind} = 0;
	$countsgood{$ind} = 0;
}

for $ind (@inds) {#for each ind
	#open infile and outfile
	$editfname = $editstem . $ind . $editex;
	unless(open(INFILE, "<", $editfname)) {
		print "Cannot open $editfname, exiting..\n";
		exit;
	}
	$neweditfname = $newdirname . '/' . $ind . $editex;
	unless(open(OUTFILE, ">", $neweditfname)) {
		print "Cannot open $neweditfname, exiting..\n";
		exit;
	}
	#read data
	while ($line1 = <INFILE>) {
		$line2 = <INFILE>;
		chomp $line1;
		chomp $line2;
		if ($line1 =~ /pair/) {#if this seq is not trimmed
			print OUTFILE "$line1\n$line2\n";
			++$countsgood{$ind};
		} else {#count a discarded read
			++$countstrim{$ind};
		}
	}
	close INFILE;
	close OUTFILE;
}
#print counts to a report file
$repfilename = $newdirname . '/' . $repfilename;
unless(open(REP, ">", $repfilename)) {
	print "Cannot open $repfilename, exiting..\n";
	exit;
}
print REP "ind\tn_discarded_trimmed_reads\tn_good_reads\n";
for $ind (@inds) {
	print REP "$ind\t$countstrim{$ind}\t$countsgood{$ind}\n";
}
close REP;

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print RUNTIMEOUT "Run took $run_s seconds.\n";
close RUNTIMEOUT;
exit;


