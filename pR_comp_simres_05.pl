#!/usr/bin/perl -w
#pR_comp_simres version 05 Copyright 2016 Andreas Hapke
#This program compares genotypes simulated with sim_fastq
#to those inferred by pyRAD version 3.0.1 based on sequence data
#produced by sim_fastq.
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

#keyword infilecolumns! means that column order in infiles matters

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
a => '',#name of *.alleles file out of pyRAD (must be defined)
s => '',#name of dir that contains files about simulated genotypes (if not defined: wd)
o => '',#name of dir for outfiles
r => ''#original read length must be the same for (f- and r-reads)
		#must be defined (positive integer); reads may not overlap
);
my %defaults = (
o => 'pR_comp_simres_out'#name of dir for outfiles
);
my $n_arg = 0;#number of command line arguments
my $flag = '';#a flag
my $val = '';#a value

#Other settings
my $repfname = 'pR_comp_simres_rep.txt';#report outfile name
my $oripopallfn = 'ori_popall.txt';#ori_popall file out of sim_fastq
my $ori_gt_fn = 'ori_gt.txt';#ori_gt file out of sim_fastq
my $loc_count_fn = 'loc_counts.txt';#name of loc_counts outfile
my $ind_gt_counts_fn = 'ind_gt_counts.txt';#name of ind_gt_counts outfile
my $perc_mis_fn = 'perc_miss_gt.txt';#name of outfile with percentages of missing genotypes
my $perc_incorr_fn = 'perc_incorr_gt.txt';#name of outfile with percentages of incorrect genotypes
my $nloc_nind_fn = 'nloc_nind.txt';#name of outfile with number of loci with at least n individuals
#variables
my $islash = '';#slash or backslash for infile path
my $oslash = "\/";#slash or backslash for outfile path
my %ilociall = ();# {ilocid}{iallid}=iallseq inferred alleles with locus and allele ids as keys
				# '-' and 'N' count as different from ACGT
my %igt = ();# {ind}{ilocID}=(iAid,iBid) inferred genotypes
my %olocoall = ();# {olocid}{oallid}=oallseq original loci and alleles
				#sub match_loci changes oallseq to reverse complement
				#where necessary to correspond to inferred alleles
my %ialloall = ();# {ilocid}{iallid}{olocid}{oallid}='f'/'r'/'b'
				#matching (dist 0 ignoring N) inferred and original alleles
				#f/r/b: original allele in f- /r- / both directions
my %ilocoloc = ();# {ilocid}=(olocid,'f'/'r') corresponding inferred and original locus ids
				#contains only correctly matching loci
				#'f','r' for orientation of oloc
my %olociloc = ();# {olocid}=ilocid contains only correctly matching loci
my %iloc_cat = ();# {ilocid}='ex'/'err'/'corr' for excluded/erroneous/correct inferred loci
my %oloc_cat = ();# {olocid}='ex'/'mis'/'pres' for excluded/missing/present original loci
my $olocid = 0;
my %ogt = ();# {ind}{olocid}=(oAid,oBid) original genotypes
my %oind = ();# {ind}=1 lasting list of original individuals
my %oallcount = ();# {oallid}=count (expected count in inferred data)
my %olocvarpos = ();# {olocid}=(array of variable positions, count starting with 0)
my %igt_cat = ();# {ind}{cat}=count cat: 'mis'/'incorr'/'corr'
				#mis: missing: genotype not observed (locus may be correct/erroneous/excluded
				#			or:locus correct, genotype matches, but
				#			not all original SNPs that are expected to occur in inferred data fully determined
				#incorr: incorrect: locus erroneous or excluded or
				#		locus correct but genotype doesn't match
				#corr: correct: locus correct, genotype matches and
				#all original SNPs that are expected to occur in inferred data fully determined
my %iloc_nind = ();# {ilocid}=nind for each iloc: number of individuals with
				#genotypes that are finally not categorized as missing (correct and incorrect genotypes)
my $success = 'y';

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
#-a must be defined
if(length $user_settings{'a'} > 0){
	unless(open(INFILE, "<", "$user_settings{'a'}")){
		print "Cannot open $user_settings{'a'}, exiting..\n";
		exit;
	}
	close INFILE;
}
else{
	print "Name of *.alleles_outfile out of pyRAD not defined (flag -a), exiting..\n";
	exit;
}
#-s if defined: Do you find the directory
if(length $user_settings{'s'} > 0){
	#remove trailing slash or backslash if any
	if($user_settings{'s'} =~ m/\//){
		$islash = "\/";
		$user_settings{'s'} =~ s/\/$//;
	}
	elsif($user_settings{'s'} =~ m/\\/){
		$islash = "\\";
		$user_settings{'s'} =~ s/\\$//;
	}
	if(length $user_settings{'s'} > 0){
		unless(-d "$user_settings{'s'}"){
			print "Cannot find infile directory $user_settings{'s'}, exiting..\n";
			exit;
		}
		#append slash or backslash to infiledir for easier construction of infilenames
		$user_settings{'s'} .= $islash;
	}	
}
#-o create if not yet existing
if (length $user_settings{'o'} > 0) {
	#remove trailing slash or backslash if any
	if($user_settings{'o'} =~ m/\//){
		$oslash = "\/";
		$user_settings{'o'} =~ s/\/$//;
	}
	elsif($user_settings{'o'} =~ m/\\/){
		$oslash = "\\";
		$user_settings{'o'} =~ s/\\$//;
	}
	if (length $user_settings{'o'} == 0){
		print "unusable outdirname, using default $defaults{'o'}\n";
		$user_settings{'o'} = $defaults{'o'};
	}
} else {
	print "no outdirname defined, using default $defaults{'o'}\n";
	$user_settings{'o'} = $defaults{'o'};
}
if (-d "$user_settings{'o'}"){
	print "Directory $user_settings{'o'} already exists.\n",
	"Please rename or delete. Exiting..\n";
	exit;
} else {
	unless(mkdir "$user_settings{'o'}"){
		print "Cannot create directory $user_settings{'o'}, exiting..\n";
		exit;
	}
}
#append slash or backslash to outdirname for easier construction of outfilenames
$user_settings{'o'} .= $oslash;

#r must be defined
unless(length $user_settings{'r'} > 0 and $user_settings{'r'} =~ /^\d+$/
and $user_settings{'r'} > 0){
	print "Original read length must be defined with flag -r (positive integer), exiting..\n";
	exit;
}
#}

####################
#Open report outfile
####################

#{
$repfname = $user_settings{'o'} . $repfname;
unless(open(OUTREP, ">", $repfname)){
	print "Cannot open $repfname, exiting..\n";
	exit;
}
#}

############
#Do the work
############

#{
#Read in *.alleles file out of pyRAD, populate %ilociall and %igt
$success = readall($user_settings{'r'},$user_settings{'a'},\%ilociall,\%igt);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}
#Read in ori_popall.txt out of sim_fastq, populate %olocoall
$oripopallfn = "$user_settings{'s'}" . $oripopallfn;
$success = read_ori_popall($oripopallfn,\%olocoall);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}
#Match inferred to original alleles
match_alleles(\%ilociall,\%olocoall,\%ialloall);

#Match loci
match_loci(\%ilociall,\%olocoall,\%ialloall,\%ilocoloc,\%olociloc,\%iloc_cat,\%oloc_cat);

#Remove missing and excluded loci from %olocoall
for $olocid (keys %oloc_cat){
	unless($oloc_cat{$olocid} eq 'pres'){
		delete $olocoall{$olocid};
	}
}
#Change sequences in %olocoall to revcomp where necessary to match corrsponding iallseqs
oallseq_rc(\%ilocoloc,\%olocoall);

#Read in ori_gt.txt, populate %ogt and %oind
$ori_gt_fn = "$user_settings{'s'}" . $ori_gt_fn;
$success = read_ori_gt($ori_gt_fn,\%ogt,\%oind);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}

#Remove unobserved and excluded genotypes from %ogt
rem_unob_ogt(\%ogt,\%oloc_cat,\%olociloc,\%igt);

#Determine expected count of each oall in inferred data (excluding excluded loci)
oall_count(\%olocoall,\%ogt,\%oallcount);

#Remove alleles that are not expected to occur in inferred data from %olocoall
$success = rem_unex_oall(\%olocoall,\%oallcount);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}
%oallcount = ();#set back, not needed any more

#Determine variable positions for each original locus
#based on those alleles that are expected to occur in inferred data
#other alleles are already removed
oloc_det_varpos(\%olocoall,\%olocvarpos);

#Compare inferred and original genotypes
comp_geno(\%igt_cat,\%iloc_nind,\%oind,\%iloc_cat,\%igt,
\%ogt,\%ilociall,\%ilocoloc,\%ialloall,\%olocvarpos);

#Produce outfile loc_counts
$loc_count_fn = $user_settings{'o'} . $loc_count_fn;
$success = loc_count_out($loc_count_fn,\%oloc_cat,\%iloc_cat);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}

#Produce outfile ind_gt_counts
$ind_gt_counts_fn = $user_settings{'o'} . $ind_gt_counts_fn;
$success = ind_gt_counts_out($ind_gt_counts_fn,\%igt_cat);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}

#Produce outfile perc_miss_gt
$perc_mis_fn = $user_settings{'o'} . $perc_mis_fn;
$success = perc_mis_out($perc_mis_fn,\%igt_cat,\%iloc_cat);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}

#Produce outfile perc_incorr_gt
$perc_incorr_fn = $user_settings{'o'} . $perc_incorr_fn;
$success = perc_incorr_out($perc_incorr_fn,\%igt_cat);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}
#Produce outfile nloc_nind
$nloc_nind_fn = $user_settings{'o'} . $nloc_nind_fn;
$success = nloc_nind_out($nloc_nind_fn,\%iloc_nind);
unless($success eq 'y'){
	print "$success\n";
	print OUTREP "$success\n";
	close OUTREP;
	exit;
}
#}

#determine runtime and print
$end_time = time();
$run_s = $end_time - $start_time;
print "Run took $run_s seconds.\n";
print OUTREP "Run took $run_s seconds.\n";
close OUTREP;

exit;

############
#Subroutines
############

#Definition of subroutine readall
#Reads in *.alleles outfile out of pyRAD
#Converts leading and trailing gaps to 'N'
#Populates datastructures owned by main:
#%ilociall = ();# {ilocid}{iallid}=iallseq inferred alleles with locus and allele ids as keys
				# '-' and 'N' count as different from ACGT
#%igt = ();# {ind}{ilocID}=(iAid,iBid) inferred genotypes
#returns 'y' if oK or error message as scalar
#expects:
#$ori_rlen original read length
#$infilename name of *.alleles outfile out of pyRAD
#ref to %ilociall
#ref to %igt
sub readall {
	#declare and initialize, _r means a reference
	my($ori_rlen, $infilename, $ilociall_r, $igt_r) = @_;
	my $seqlen = 0;
	my $success = 'y';
	my $inline1 = '';
	my $inline2 = '';
	my $fread = '';
	my $rread = '';
	my $ind = '';#ind id
	my $ilocid = 1;#inferred locus id
	my $iallid = 0;#inferred allele id
	my $new_iallid = 1;#new inferred allele id
	my $iAid = 0;#first allele id of genotype
	my $iBid = 0;#2nd allele id of genotype
	my $iallseq = '';#inferred allele sequence
	my $pos = 0;
	my $tempstring1 = '';
	my @temparr1 = ();
	my @temparr2 = ();
	my %temphash1 = ();
	
	print "Reading $infilename...\n";	
	#open infile and read in
	unless(open(INFILE, "<", $infilename)){
		$success = "Cannot open $infilename.";
		return $success;
	}
	while($inline1 = <INFILE>){
		$inline1 =~ s/(\n|\r\n?)//;#should remove windows and linux line endings
		if($inline1 =~ /^>/){#Line contains a header and sequence
			$inline2 = <INFILE>;#get 2nd sequence of this ind
			$inline2 =~ s/(\n|\r\n?)//;#should remove windows and linux line endings
			unless($inline2 =~ /^>/){
				$success = "implausible line in $infilename:\n$inline2";
				return $success;
			}
			#get first allele
			$inline1 =~ s/^>//;#remove leading char of header
			@temparr1 = split(/\s+/,$inline1);
			@temparr2 = split(/_/,$temparr1[0]);
			$ind = $temparr2[0];
			$iallseq = $temparr1[1];
			#Check length of iallseq
			if((length $iallseq) < (2 * $ori_rlen + 4)){#characters are missing in sequence
				#This happens in loci with few individuals where the first position is always missing.
				#Append N to beginning of sequence until it has sufficient length.
				while((length $iallseq) < (2 * $ori_rlen + 4)){
					$iallseq = 'N' . $iallseq;
				}
			}
			unless((length $iallseq) >= (2 * $ori_rlen + 4)){
				$success = "Implausible sequence length in $infilename. Sequence:\n$iallseq\n";
				return $success;
			}
			#Convert leading and trailing gaps to 'N'
			for ($pos = 0; $pos < length $iallseq; ++$pos){
				if(substr($iallseq,$pos,1) eq '-'){
					substr($iallseq,$pos,1) = 'N';				
				} else {
					last;				
				}
			}
			for ($pos = (length $iallseq) - 1; $pos >= 0; --$pos){
				if(substr($iallseq,$pos,1) eq '-'){
					substr($iallseq,$pos,1) = 'N';				
				} else {
					last;				
				}
			}
			#Get sequence without spacer and additional bases adjacent to spacer
			$fread = substr($iallseq,0,$ori_rlen);
			$rread = substr($iallseq,((length $iallseq) - $ori_rlen),$ori_rlen);
			$iallseq = $fread . $rread;
			if(defined $temphash1{$iallseq}){#already known allele sequence
				$iAid = $temphash1{$iallseq};
			} else {#new allele sequence
				$iAid = $new_iallid;
				$temphash1{$iallseq} = $iAid;
				++$new_iallid;
			}
			#get 2nd allele
			$inline2 =~ s/^>//;#remove leading char of header
			@temparr1 = split(/\s+/,$inline2);
			@temparr2 = split(/_/,$temparr1[0]);
			unless($ind eq $temparr2[0]){
				$success = "implausible line in $infilename:\n$inline2";
				return $success;
			}
			$iallseq = $temparr1[1];
			#Check length of iallseq
			if((length $iallseq) < (2 * $ori_rlen + 4)){#characters are missing in sequence
				#This happens in loci with few individuals where the first position is always missing.
				#Append N to beginning of sequence until it has sufficient length.
				while((length $iallseq) < (2 * $ori_rlen + 4)){
					$iallseq = 'N' . $iallseq;
				}
			}
			unless(length $iallseq >= (2 * $ori_rlen + 4)){
				$success = "Implausible sequence length in $infilename. Sequence:\n$iallseq\n";
				return $success;
			}
			#Convert leading and trailing gaps to 'N'
			for ($pos = 0; $pos < length $iallseq; ++$pos){
				if(substr($iallseq,$pos,1) eq '-'){
					substr($iallseq,$pos,1) = 'N';				
				} else {
					last;				
				}
			}
			for ($pos = (length $iallseq) - 1; $pos >= 0; --$pos){
				if(substr($iallseq,$pos,1) eq '-'){
					substr($iallseq,$pos,1) = 'N';				
				} else {
					last;				
				}
			}
			#Get sequence without spacer and additional bases adjacent to spacer
			$fread = substr($iallseq,0,$ori_rlen);
			$rread = substr($iallseq,((length $iallseq) - $ori_rlen),$ori_rlen);
			$iallseq = $fread . $rread;
			if(defined $temphash1{$iallseq}){#already known allele sequence
				$iBid = $temphash1{$iallseq};
			} else {#new allele sequence
				$iBid = $new_iallid;
				$temphash1{$iallseq} = $iBid;
				++$new_iallid;
			}
			#store genotype in %igt owned by main
			@{$$igt_r{$ind}{$ilocid}} = ($iAid,$iBid);			
		}
		elsif($inline1 =~ m{^\/\/}){#end of locus reached
			#store alleles collected in %temphash1 in %ilociall owned by main
			for $iallseq(keys %temphash1){
				$iallid = $temphash1{$iallseq};
				$$ilociall_r{$ilocid}{$iallid}=$iallseq;
			}
			%temphash1 = ();#set back
			++$ilocid;#increment for next locus
		}
		else{
			$success = "implausible line in $infilename:\n$inline1";
			return $success;
		}
	}
	close INFILE;
	return $success;
}

#Definition of subroutine read_ori_popall
#Reads in ori_popall.txt out of sim_dat
#Populates %olocoall owned by main:
# {olocid}{oallid}=oallseq original loci and alleles
#returns 'y' if oK or error message as scalar
#expects:
#infilename
#ref to %olocoall
sub read_ori_popall{
	#declare and initialize, _r means a reference
	my($infilename, $olocoall_r) = @_;
	my $inline = '';
	my @temparr1 = ();
	my $success = 'y';
	
	print "Reading $infilename...\n";	
	unless(open(INFILE, "<", $infilename)){
		$success = "Cannot open $infilename";
		return $success;
	}
	$inline = <INFILE>;#skip header line
	while($inline = <INFILE>){
		$inline =~ s/(\n|\r\n?)//;#should remove windows and linux line endings
		@temparr1 = split(/\t/,$inline);
		#infilecolumns!
		$$olocoall_r{$temparr1[0]}{$temparr1[1]} = $temparr1[3];
	}
	close INFILE;
	return $success;
}

#Definition of subroutine match_alleles
#Matches inferred to original alleles
#populates %ialloall owned by main:
# {ilocid}{iallid}{olocid}{oallid}='f'/'r'/'b'
#expects:
#ref to %ilociall {ilocid}{iallid}=iallseq
#ref to %olocoall {olocid}{oallid}=oallseq
#ref to %ialloall see above
sub match_alleles{
	#declare and initialize, _r means a reference
	my($ilociall_r,$olocoall_r,$ialloall_r) = @_;
	my %oall = ();# {oallseq}=(olocid,oallid)
	my $ilocid = 0;
	my $iallid = 0;
	my $iallseq = '';
	my $riallseq = '';#
	my $olocid = 0;
	my $oallid = 0;
	my $folocid = 0;#f-matching olocid
	my $rolocid = 0;#r-matching olocid
	my $foallid = 0;#f-matching oallid
	my $roallid = 0;#r-matching oallid
	my $oallseq = '';
	my $nN = 0;#number of N in a sequence
	my $Nmaxperm = 3;#maximum number of N in inferred sequence for permutation method
	my @Nallperm = ('600','2400','10000','40000');
	my @iallper = ();#all possible permutations of a sequence containing N (replacing N with ACGT)
	my $noriall = 0;#number of original alleles	
	my $fmatch = 0;
	my $rmatch = 0;
	my $niloc = 0;
	my $nilocdone = 0;
	my $i = 0;
	
	print "Matching alleles...\n";	
	#Just for reporting progress
	$niloc = keys %{$ilociall_r};
	#populate %oall
	for $olocid (keys %{$olocoall_r}){
		for $oallid (keys %{$$olocoall_r{$olocid}}){
			$oallseq = $$olocoall_r{$olocid}{$oallid};
			@{$oall{$oallseq}} = ($olocid,$oallid);
		}
	}
	#adjust $nmaxperm according to number of original alleles
	$noriall = keys %oall;	
	for $i (@Nallperm){
		if($noriall >= $i){
			++$Nmaxperm;	
		} else {
			last;		
		}	
	}
	#match inferred alleles to original alleles
	for $ilocid (keys %{$ilociall_r}){
		for $iallid (keys %{$$ilociall_r{$ilocid}}){
			$iallseq = $$ilociall_r{$ilocid}{$iallid};
			if($iallseq =~ m/N/){
				$nN = $iallseq =~ tr/N//;#count N in sequence
				if($nN <= $Nmaxperm){#sequence has less than $Nmaxperm N
					#build all possible permutations of sequence and lookup in %oall
					@iallper = nper($iallseq);
					#search in f-orientation
					for $iallseq (@iallper){
						if(defined $oall{$iallseq}){
							($olocid,$oallid) = @{$oall{$iallseq}};
							$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid} = 'f';
						}
					}
					#search in r-orientation
					for $iallseq (@iallper){
						$riallseq = reverse $iallseq;
						$riallseq =~ tr/ACGTacgt/TGCAtgca/;
						if(defined $oall{$riallseq}){
							($olocid,$oallid) = @{$oall{$riallseq}};
							#check if any permutation of the same iallseq matched in f-direction (improbable...)
							#many 'if defined' but necessary to avoid autovivification of hash keys
							if(defined $$ialloall_r{$ilocid} and defined $$ialloall_r{$ilocid}{$iallid}
							and defined $$ialloall_r{$ilocid}{$iallid}{$olocid} and defined
							$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid}){
								$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid} = 'b';
							} else {
								$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid} = 'r';
							}
						}
					}
				}
				else {#sequence has 5 or more N, compare to each oallseq
					$riallseq = reverse $iallseq;
					$riallseq =~ tr/ACGTacgt/TGCAtgca/;
					for $oallseq (keys %oall){
						#check forward
						if(distpwdN($iallseq,$oallseq) == 0){
							$fmatch = 1;
							($olocid,$oallid) = @{$oall{$oallseq}};
						}
						#check reverse complement
						if(distpwdN($riallseq,$oallseq) == 0){
							$rmatch = 1;
							($olocid,$oallid) = @{$oall{$oallseq}};
						}
						if($fmatch == 1 and $rmatch == 0){#match in f-orientation only
							$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid} = 'f';
							$fmatch = 0;
							$olocid = 0;
							$oallid = 0;
						}
						elsif($fmatch == 0 and $rmatch == 1){#match in r-orientation only
							$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid} = 'r';
							$rmatch = 0;
							$olocid = 0;
							$oallid = 0;
						}
						elsif($fmatch == 1 and $rmatch == 1){#match in both orientations
							$$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid} = 'b';
							$fmatch = 0;
							$rmatch = 0;
							$olocid = 0;
							$oallid = 0;							
						}
					}
				}
			}
			else {#$iallseq has no N
				#lookup in f-orientation
				if(defined $oall{$iallseq}){
					$fmatch = 1;
					($folocid,$foallid) = @{$oall{$iallseq}};
				}
				#lookup in r-orientation
				$riallseq = reverse $iallseq;
				$riallseq =~ tr/ACGTacgt/TGCAtgca/;
				if(defined $oall{$riallseq}){
					$rmatch = 1;
					($rolocid,$roallid) = @{$oall{$riallseq}};
				}
				if($fmatch == 1 and $rmatch == 0){#match in f-orientation only
					$$ialloall_r{$ilocid}{$iallid}{$folocid}{$foallid} = 'f';
					$fmatch = 0;
					$folocid = 0;
					$foallid = 0;
				}
				elsif($fmatch == 0 and $rmatch == 1){#match in r-orientation only
					$$ialloall_r{$ilocid}{$iallid}{$rolocid}{$roallid} = 'r';
					$rmatch = 0;
					$rolocid = 0;
					$roallid = 0;
				}
				elsif($fmatch == 1 and $rmatch == 1){#match in both orientations
					if($foallid == $roallid){#iallseq matches same oallseq as f and rc (improbable)
											#allele ids are unique, sufficient to compare allele id
						$$ialloall_r{$ilocid}{$iallid}{$folocid}{$foallid} = 'b';
					} else {
						$$ialloall_r{$ilocid}{$iallid}{$folocid}{$foallid} = 'f';
						$$ialloall_r{$ilocid}{$iallid}{$rolocid}{$roallid} = 'r';
					}
					$fmatch = 0;
					$rmatch = 0;
					$folocid = 0;
					$foallid = 0;
					$rolocid = 0;
					$roallid = 0;
				}
			}
		}
		++$nilocdone;
		print "matched alleles for $nilocdone of $niloc loci.\n";
	}	
}

#Definition of subroutine nper
#creates all possible permutations of a DNA-sequence containing at least 1 N.
#Permutations: N replaced by A/C/G/T
#expects:
#$seq	the sequence
#returns @per: array containing all permutations
sub nper {
	#declare and initialize
	my ($seq) = @_;
	my $seql = 0;
	my $Nstring = '';#consists of $seql 'N'
	my $Nmatch = '';#contains \0 at all positions where $seq has N
					#defined values at other positions
	my @per = ();#will hold the permutations
	my $nseq = 0;#number of sequences currently in @per
	my $newper = '';#one sequence that is currently permuted	
	my $pos = 0;#current sequence position
	my $i  = 0;
	
	$seql = length $seq;
	$Nstring = "N" x $seql;#create string of seqlength x N
	$Nmatch = $seq ^ $Nstring;#see above
	$seq =~ tr/N/A/;#transliterate all N in seq by A
	$per[0] = $seq;#put into @per as first permutation
	#loop through positions
	for ($pos = 0; $pos < $seql; ++$pos) {
		#if this is an N-position
		if (substr($Nmatch,$pos,1) =~ /\0/) {
			#determine number of seqs currently in @per
			$nseq = @per;
			#loop through these seqs
			for ($i = 0; $i < $nseq; ++$i) {
				$newper = $per[$i];#copy current seq
				substr ($newper,$pos,1) = 'C';#replace current pos with C
				push @per, $newper;#put the new permutation into @per
				substr ($newper,$pos,1) = 'G';#replace current pos with G
				push @per, $newper;#put the new permutation into @per
				substr ($newper,$pos,1) = 'T';#replace current pos with T
				push @per, $newper;#put the new permutation into @per
			}
		}
	}
	return @per;
}

#definition of subroutine distpwdN version 06
#compares 2 DNA sequences (must be of equal length!)
#treats N as equal with -ACGT
#knows chars: : -ACGTN
#expects
#$seq1, $seq2
#returns 0 if sequences match and 1 if they are different
sub distpwdN {
	#declare and initialize
	my ($seq1,$seq2) = @_;
	my $seql = 0;
	my $dist = 0;#the result: 0 for match, 1 for difference
	my $Nstring = '';#a string of N of length seqlength
	my $NAseq1 = '';#seq1 with -CGT replaced by AAAA
	my $NAseq2 = '';#seq2 with -CGT replaced by AAAA
	my $changepos = '';#string that contains defined values at positions
						#where only one of the 2 seqs has N, other pos: \0
	my$pos = 0;#a position;
	
	$seql = length $seq1;
	if($seq1 =~ /N/ or $seq2 =~ /N/){
		#replace -ACGT with N at all positions that have N in any sequence
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
	}
	if($seq1 eq $seq2){
		$dist = 0;
	} else {
		$dist = 1;
	}
	return $dist;
}

#Definition of subroutine match_loci
#Finds correctly matching loci, excluded, erroneous, and missing loci
#Correct match: one iloc matches one oloc
#with only f and/or b, r and/or b, or b oriented matches
#expects
#ref to %ilociall {ilocid}{iallid}=iallseq
#ref to %olocoall {olocid}{oallid}=oallseq
#ref to %ialloall {ilocid}{iallid}{olocid}{oallid}='f'/'r'/'b'
#ref to %ilocoloc {ilocid}=(olocid,'f'/'r')
#ref to %olociloc {olocid}=ilocid
#ref to %iloc_cat {ilocid}='ex'/'err'/'corr'
#ref to %oloc_cat {olocid}='ex'/'mis'/'pres'
#populates %ilocoloc, %olociloc, %iloc_cat, %oloc_cat
sub match_loci{
	#declare and initialize, _r means a reference
	my($ilociall_r,$olocoall_r,$ialloall_r,$ilocoloc_r,$olociloc_r,$iloc_cat_r,$oloc_cat_r) = @_;
	my $ilocid = 0;
	my $iallid = 0;
	my $olocid = 0;
	my $oallid = 0;
	my $orient = '';# f/r/b orientation of match
	my %orients = ();# count match orientations for one iloc
	my %match_oloc = ();#matching olocs
	my @match_iloc = ();#matching ilocs
	my %olociloc = ();# {olocid}{ilocid} = 1
	
	print "Matching loci...";	
	#Match ilocs to olocs
	for $ilocid (keys %{$ilociall_r}){
		if(defined $$ialloall_r{$ilocid}){#if this iloc matches an oloc
			for $iallid (keys %{$$ialloall_r{$ilocid}}){
				for $olocid (keys %{$$ialloall_r{$ilocid}{$iallid}}){
					$match_oloc{$olocid} = 1;
					for $oallid (keys %{$$ialloall_r{$ilocid}{$iallid}{$olocid}}){
						$orient = $$ialloall_r{$ilocid}{$iallid}{$olocid}{$oallid};
						++$orients{$orient};
					}
				}
			}
			if(keys %match_oloc > 1){#iloc matches > 1 oloc
				$$iloc_cat_r{$ilocid} = 'ex';#exclude iloc
				for $olocid (keys %match_oloc){
					$$oloc_cat_r{$olocid} = 'ex';#exclude all matching olocs
				}
			} else {#iloc matches one oloc
				$olocid = (keys %match_oloc)[0];
				if(defined $orients{'f'} and defined $orients{'r'}){#f- and r-match
					#break					
					#print "f- and r-match\n";
					#endbreak
					$$iloc_cat_r{$ilocid} = 'ex';#exclude iloc
					$$oloc_cat_r{$olocid} = 'ex';#exclude oloc
				}
				elsif(defined $orients{'r'}){#r-match or r- and b-match
					$$iloc_cat_r{$ilocid} = 'corr';
					$$oloc_cat_r{$olocid} = 'pres';
					@{$$ilocoloc_r{$ilocid}} = ($olocid,'r');
				} else {#f-match, f- and b- match or b-match
					$$iloc_cat_r{$ilocid} = 'corr';
					$$oloc_cat_r{$olocid} = 'pres';
					@{$$ilocoloc_r{$ilocid}} = ($olocid,'f');
				}
			}
			%match_oloc = ();#set back
			%orients = ();#set back
		}
		else{#iloc matches no oloc
			$$iloc_cat_r{$ilocid} = 'err';#iloc is erroneous
		}
	}
	#populate %olociloc
	for $ilocid (keys %{$ialloall_r}){
		for $iallid (keys %{$$ialloall_r{$ilocid}}){
			for $olocid (keys %{$$ialloall_r{$ilocid}{$iallid}}){
				$olociloc{$olocid}{$ilocid} = 1;
			}
		}
	}
	#match olocs to ilocs
	for $olocid (keys %{$olocoall_r}){
		if(defined $olociloc{$olocid}){#oloc matches an iloc
			for $ilocid (keys %{$olociloc{$olocid}}){
				push @match_iloc, $ilocid;
			}
			if(@match_iloc > 1){#more than one iloc matches
				$$oloc_cat_r{$olocid} = 'ex';#exclude oloc
				for $ilocid (@match_iloc){
					$$iloc_cat_r{$ilocid} = 'ex';#exclude all matching ilocs
				}
			}
			@match_iloc = ();#set back
		} else {#oloc matches no iloc
			$$oloc_cat_r{$olocid} = 'mis';#oloc is missing
		}
	}
	#remove ilocs that are now flagged 'ex' from %ilocoloc owned by main
	for $ilocid (keys %{$iloc_cat_r}){
		if($$iloc_cat_r{$ilocid} eq 'ex'){
			if(defined $$ilocoloc_r{$ilocid}){
				delete $$ilocoloc_r{$ilocid};
			}
		}
	}
	#populate %olociloc owned by main
	for $ilocid (keys %{$ilocoloc_r}){
		$olocid = $$ilocoloc_r{$ilocid}[0];
		$$olociloc_r{$olocid} = $ilocid;
	}
}

#Definition of sub oallseq_rc
#Changes sequences in %olocoall owned by main to reverse complement
#where necessary to match corresponding iallseqs
#expects
#ref to %ilocoloc {ilocid}=(olocid,'f'/'r') r: revcomp of oloc matches iloc
#ref to %olocoall {olocid}{oallid}=oallseq
sub oallseq_rc{
	#declare and initialize, _r means a reference
	my($ilocoloc_r,$olocoall_r) = @_;
	my $ilocid = 0;
	my $olocid = 0;
	my $orient = '';
	my $oallid = 0;
	my $oallseq = '';
	
	for $ilocid (keys %{$ilocoloc_r}){
		($olocid,$orient) = @{$$ilocoloc_r{$ilocid}};
		if($orient eq 'r'){
			for $oallid (keys %{$$olocoall_r{$olocid}}){
				$oallseq = $$olocoall_r{$olocid}{$oallid};
				$oallseq = reverse $oallseq;
				$oallseq =~ tr/ACGTacgt/TGCAtgca/;
				$$olocoall_r{$olocid}{$oallid} = $oallseq;
			}
		}
	}
}

#Definition of sub read_ori_gt
#Reads in ori_gt.txt out of sim_dat
#Populates data structures owned by main %ogt and %oind
#returns 'y' if successful or error message as scalar
#expects
#$infilename
#ref to %ogt {ind}{olocid}=(oAid,oBid)
#ref to %oind {ind}=1
sub read_ori_gt{
	#declare and initialize, _r means a reference
	my($infilename,$ogt_r,$oind_r) = @_;
	my $success = 'y';
	my $inline = '';
	my @temparr1 = ();
	my $ind = '';
	
	print "Reading $infilename...\n";	
	unless(open(INFILE, "<", $infilename)){
		$success = "Cannot open $infilename.";
		return $success;
	}
	$inline = <INFILE>;#skip header line
	while($inline = <INFILE>){
		$inline =~ s/(\n|\r\n?)//;#should remove windows and linux line endings
		@temparr1 = split(/\t/,$inline);
		#infilecolumns!
		@{$$ogt_r{$temparr1[0]}{$temparr1[1]}} = ($temparr1[2],$temparr1[3]);
	}
	close INFILE;
	for $ind (keys %{$ogt_r}){
		$$oind_r{$ind} = 1;
	}
	return $success;
}

#Definition of sub rem_unob_ogt
#Removes unobserved and excluded genotypes from %ogt owned by main
#expects
#ref to %ogt {ind}{olocid}=(oAid,oBid) original genotypes
#ref to %oloc_cat {olocid}='ex'/'mis'/'pres' for excluded/missing/present original loci
#ref to %olociloc {olocid}=ilocid contains only correctly matching loci
#ref to %igt {ind}{ilocID}=(iAid,iBid) inferred genotypes
sub rem_unob_ogt{
	#declare and initialize, _r means a reference
	my($ogt_r,$oloc_cat_r,$olociloc_r,$igt_r) = @_;
	my $ind = '';
	my $olocid = 0;
	my $ilocid = 0;
	
	for $ind (keys %{$ogt_r}){
		for $olocid (keys %{$$ogt_r{$ind}}){
			if($$oloc_cat_r{$olocid} eq 'pres'){#locus is expected in inferred data
				$ilocid = $$olociloc_r{$olocid};
				unless(defined $$igt_r{$ind} and defined $$igt_r{$ind}{$ilocid}){
					delete $$ogt_r{$ind}{$olocid};
				}
			} else {#locus is missing in inferred data or excluded
				delete $$ogt_r{$ind}{$olocid};
			}
		}
	}
	#Delete empty inds from %ogt owned by main if any
	for $ind (keys %{$ogt_r}){
		if(keys %{$$ogt_r{$ind}} == 0){
			delete $$ogt_r{$ind};
		}
	}
}

#Definition of sub oall_count
#Determines expected count of oalls in inferred data
#excluding excluded loci, based on original genotypes
#now, these contain only genotypes that have corresponding genotypes in inferred data
#expects
#ref to %olocoall {olocid}{oallid}=oallseq
#ref to %ogt {ind}{olocid}=(oAid,oBid)
#ref to %oallcount {oallid}=count
sub oall_count{
	#declare and initialize, _r means a reference
	my($olocoall_r,$ogt_r,$oallcount_r) = @_;
	my $olocid = 0;
	my $oallid = 0;
	my $ind = '';
	my $oAid = 0;#first allele of genotype
	my $oBid = 0;#second id of genotpye
	
	#initialize allele counts with 0
	for $olocid (keys %{$olocoall_r}){
		for $oallid (keys %{$$olocoall_r{$olocid}}){
			$$oallcount_r{$oallid} = 0;
		}
	}
	#count
	for $ind (keys %{$ogt_r}){
		for $olocid (keys %{$$ogt_r{$ind}}){
			($oAid,$oBid) = @{$$ogt_r{$ind}{$olocid}};
			++$$oallcount_r{$oAid};
			++$$oallcount_r{$oBid};
		}
	}
}

#Definition of sub rem_unex_oall
#Removes alleles from %olocoall owned by main
#that are not expected to occur in inferred data
#returns 'y' or error message as scalar
#expects:
#ref to %olocoall {olocid}{oallid}=oallseq
#ref to %oallcount {oallid}=count
sub rem_unex_oall{
	#declare and initialize, _r means a reference
	my($olocoall_r,$oallcount_r) = @_;
	my $olocid = 0;
	my $oallid = 0;
	my $success = 'y';
	
	for $olocid (keys %{$olocoall_r}){
		for $oallid (keys %{$$olocoall_r{$olocid}}){
			if($$oallcount_r{$oallid} == 0){
				delete $$olocoall_r{$olocid}{$oallid};
			}
		}
	}
	#report error if any loci in %olocoall owned by main are empty now
	for $olocid (keys %{$olocoall_r}){
		if(keys %{$$olocoall_r{$olocid}} == 0){
			$success = "rem_unex_oall: found empty locus in %olocoall.";
			return $success;
		}
	}
	return $success;
}

#Definition of subroutine oloc_det_varpos
#Determines variable positions for each locus in %olocoall
#populates %olocvarpos owned by main
#expects
#ref to %olocoall {olocid}{oallid}=oallseq
#ref to %olocvarpos {olocid}=(array of variable positions, count starting with 0)
sub oloc_det_varpos{
	#declare and initialize, _r means a reference
	my($olocoall_r,$olocvarpos_r) = @_;
	my $olocid = 0;
	my $oallid = 0;
	my $oallseq = '';
	my @seqvar = ();#one allele sequence split into array
	my @oallseqvar = ();#2d array of sequences of one locus
					#d1: seq, d2: pos
	my @varpos = ();
	
	for $olocid (keys %{$olocoall_r}){
		for $oallid (keys %{$$olocoall_r{$olocid}}){
			$oallseq = $$olocoall_r{$olocid}{$oallid};
			@seqvar = split('',$oallseq);
			push @oallseqvar, [@seqvar];
		}
		@varpos = varpos(\@oallseqvar);
		@{$$olocvarpos_r{$olocid}} = @varpos;
		@oallseqvar = ();#set back
	}
}

#Definition of subroutine varpos
#expects: ref to sequence alignment as 2d array:
#d1: rows, d2: cols
#returns: an array of variable positions
#if no positions are variable: empty
#allowed characters: ACGTN, only capital, N is a 5th char
#all sequences must have same length
sub varpos {
	#declare and initialize
	my ($seqmat_r) = @_;
	my @seqmat = @{$seqmat_r};
	my @varpos = ();#the variable positions
	my $row = 0;
	my $col = 0;
	my $char1 = '';#a character in the first sequence
	my %count = ();#a hash for counting characters
	
	#loop through cols
	for ($col = 0; $col < @{$seqmat[0]}; ++$col) {
		$char1 = ${$seqmat[0]}[$col];#character of 1st sequence in this col
		#count how many sequences have this character:
		#loop through rows
		for ($row = 0; $row < @seqmat; ++$row) {
			++$count{${$seqmat[$row]}[$col]};
		}
		#if not all have this character:
		if ($count{$char1} < @seqmat) {
			#add this col to @varpos
			push @varpos, $col;
		}
		%count = ();#set counthash back
	}
	return @varpos;
}

#Definition of sub comp_geno
#compares inferred to original genotypes
#populates %igt_cat and %iloc_nind owned by main
#Categories:
#mis: missing: genotype not observed or
#             genotype observed and matching but not all original SNPs fully determined
#			  original SNPs include only those expected to occur in inferred dataset
#incorr: incorrect: locus is excluded or erroneous or
#					locus is correct but genotype doesn't match
#corr: locus is correct, genotype matches, and all original SNPs are fully determined (see above)
#expects references to:
#%igt_cat {ind}{cat}=count cat:'mis'/'incorr'/'corr' see above
#%iloc_nind {ilocid}=nind for each iloc: number of individuals with
					#genotypes that are finally not categorized as missing (correct and incorrect genotypes)
#%oind {ind}=1 complete list of all original inds
#%iloc_cat {ilocid}='ex'/'err'/'corr'
#%igt {ind}{ilocID}=(iAid,iBid) inferred genotypes
#%ogt {ind}{olocid}=(oAid,oBid) original genotypes
#%ilociall {ilocid}{iallid}=iallseq
#%ilocoloc {ilocid}=(olocid,'f'/'r') contains only correctly matching loci
#%ialloall {ilocid}{iallid}{olocid}{oallid}='f'/'r'/'b' matching alleles
#%olocvarpos {olocid}=(array of variable positions, count starting with 0)
sub comp_geno{
	#declare and initialize, _r means a reference
	my($igt_cat_r,$iloc_nind_r,$oind_r,$iloc_cat_r,$igt_r,$ogt_r,
	$ilociall_r,$ilocoloc_r,$ialloall_r,$olocvarpos_r) = @_;
	my $ind = '';#individual id
	my $ilocid = 0;
	my $iallid = 0;
	my $olocid = 0;
	my $oallid = 0;
	my $iAid = 0;#inferred genotype: id of allele A
	my $iBid = 0;#inferred genotype: id of allele B
	my $iAseq = '';
	my $iBseq = '';
	my %oallid_match_iA = ();# {oallid}=1, ids of oalleles that match inferred allele A
	my %oallid_match_iB = ();# {oallid}=1, ids of oalleles that match inferred allele B
	my $oAid = 0;#original genotype: id of allele A
	my $oBid = 0;#original genotype: id of allele B
	my $pos = 0;
	my $undet = 0;#undetermined == 1
	
	print "Comparing genotypes...\n";	
	#initialize all counts in %igt_cat owned by main with 0
	for $ind (keys %{$oind_r}){
		$$igt_cat_r{$ind}{'mis'} = 0;
		$$igt_cat_r{$ind}{'incorr'} = 0;
		$$igt_cat_r{$ind}{'corr'} = 0;		
	}
	#initialize all counts in %iloc_nind owned by main with 0
	for $ilocid (keys %{$iloc_cat_r}){
		$$iloc_nind_r{$ilocid} = 0;
	}
	#compare genotypes
	for $ind (keys %{$igt_cat_r}){
		for $ilocid (keys %{$iloc_cat_r}){
			if(defined $$igt_r{$ind}
				and defined $$igt_r{$ind}{$ilocid}){#genotype observed
				if($$iloc_cat_r{$ilocid} eq 'corr'){#locus is correct
					($iAid,$iBid) = @{$$igt_r{$ind}{$ilocid}};#get inferred allele ids
					#collect ids of original alleles that match iA and iB
					for $olocid (keys %{$$ialloall_r{$ilocid}{$iAid}}){
						for $oallid (keys %{$$ialloall_r{$ilocid}{$iAid}{$olocid}}){
							$oallid_match_iA{$oallid} = 1;
						}
					}
					for $olocid (keys %{$$ialloall_r{$ilocid}{$iBid}}){
						for $oallid (keys %{$$ialloall_r{$ilocid}{$iBid}{$olocid}}){
							$oallid_match_iB{$oallid} = 1;
						}
					}
					#get original genotype
					$olocid = $$ilocoloc_r{$ilocid}[0];
					($oAid,$oBid) = @{$$ogt_r{$ind}{$olocid}};
					if((defined $oallid_match_iA{$oAid}
						and defined $oallid_match_iB{$oBid})
						or (defined $oallid_match_iA{$oBid}
						and defined $oallid_match_iB{$oAid})){#genotypes match
						$iAseq = $$ilociall_r{$ilocid}{$iAid};
						$iBseq = $$ilociall_r{$ilocid}{$iBid};
						if($iAseq =~ /N/ or $iBseq =~ /N/){#any seq not fully determined?
							#check if all original SNPs are fully determined
							for $pos (@{$$olocvarpos_r{$olocid}}){
								if(substr($iAseq,$pos,1) eq 'N'
									or substr($iBseq,$pos,1) eq 'N'){#a SNP pos not fully determined
									$undet = 1;
									last;
								}
							}	
							if($undet > 0){#at least one SNP position undetermined
								++$$igt_cat_r{$ind}{'mis'};#count missing genotype
								$undet = 0;#set back
							} else {#all SNP positions fully determined
								++$$igt_cat_r{$ind}{'corr'};#count correct genotype
								++$$iloc_nind_r{$ilocid};#count individual with data for this locus
							}
						} else {#bot iseqs fully determined
							++$$igt_cat_r{$ind}{'corr'};#count correct genotype
							++$$iloc_nind_r{$ilocid};#count individual with data for this locus
						}
					} else {#genotypes don't match
						++$$igt_cat_r{$ind}{'incorr'};#count incorrect genotype
						++$$iloc_nind_r{$ilocid};#count individual with data for this locus
					}
					%oallid_match_iA = ();#set back
					%oallid_match_iB = ();#set back
				} else {#locus is erroneous or excluded
					++$$igt_cat_r{$ind}{'incorr'};#count incorrect genotype
					++$$iloc_nind_r{$ilocid};#count individual with data for this locus
				}
			} else {#genotype not observed
				++$$igt_cat_r{$ind}{'mis'};#count missing genotype
			}
		}
	}
}

#Definition of sub loc_count_out
#Determines counts of original and inferred locus of different categories
#Produces outfile
#expects
#$outfn outfilename
#ref to %oloc_cat {olocid}='ex'/'mis'/'pres'
#ref to %iloc_cat {ilocid}='ex'/'err'/'corr'
#returns 'y' or error message as scalar
sub loc_count_out{
	#declare and initialize, _r means a reference
	my($outfn,$oloc_cat_r,$iloc_cat_r) = @_;
	my $noloc = 0;#total number of original loci
	my $niloc = 0;#total number of inferred loci
	my %oloc_count = ();#counts of categories
	my %iloc_count = ();#counts of categories
	my $cat = '';
	my $olocid = 0;
	my $ilocid = 0;
	my $success = 'y';

	$noloc = keys %{$oloc_cat_r};
	$niloc = keys %{$iloc_cat_r};
	
	#initialize counts
	$oloc_count{'ex'} = 0;
	$oloc_count{'mis'} = 0;
	$oloc_count{'pres'} = 0;
	$iloc_count{'ex'} = 0;
	$iloc_count{'err'} = 0;
	$iloc_count{'corr'} = 0;
	
	#count original loci of different categories
	for $olocid (keys %{$oloc_cat_r}){
		$cat = $$oloc_cat_r{$olocid};
		++$oloc_count{$cat};
	}
	#count inferred loci of different categories
	for $ilocid (keys %{$iloc_cat_r}){
		$cat = $$iloc_cat_r{$ilocid};
		++$iloc_count{$cat};
	}
	unless(open(OUTFILE, ">", $outfn)){
		$success = "Cannot open outfile $outfn.";
		return $success;
	}
	print OUTFILE "original loci total\t$noloc\n",
	"found original loci\t$oloc_count{'pres'}\n",
	"excluded original loci\t$oloc_count{'ex'}\n",
	"missing original loci\t$oloc_count{'mis'}\n",
	"inferred loci total\t$niloc\n",
	"correctly inferred loci\t$iloc_count{'corr'}\n",
	"excluded inferred loci\t$iloc_count{'ex'}\n",
	"erroneous inferred loci\t$iloc_count{'err'}\n";
	close OUTFILE;	
	return $success;
}

#Definition of subroutine ind_gt_counts_out
#Produces outfile with genotype counts of different categories for each ind
#Categories: mis (missing), incorr (incorrect), corr (correct)
#returns 'y' or error message as scalar
#expects
#$outfn outfilename
#ref to %igt_cat {ind}{cat}=count cat: 'mis'/'incorr'/'corr'
sub ind_gt_counts_out{
	#Declare and initialize, _r means a reference
	my($outfn,$igt_cat_r) = @_;
	my $ind = '';
	my $success = 'y';
	
	unless(open(OUTFILE, ">", $outfn)){
		$success = "Cannot open outfile $outfn.";
		return $success;
	}
	#print header line
	print OUTFILE "ind\tcorrect\tincorrect\tmissing\n";
	#print data
	for $ind (sort keys %{$igt_cat_r}){
		print OUTFILE "$ind\t",
		"$$igt_cat_r{$ind}{'corr'}\t",
		"$$igt_cat_r{$ind}{'incorr'}\t",
		"$$igt_cat_r{$ind}{'mis'}\n";		
	}
	close OUTFILE;
	return $success;
}

#Definition of subroutine perc_mis_out
#Determines percentage of missing genotypes for each ind
#based on total number of loci in inferred dataset
#Determines avg, min, max across inds
#Produces outfile with these values
#returns 'y' or error message as scalar
#expects
#$outfn outfilename
#ref to %igt_cat {ind}{cat}=count cat: 'mis'/'incorr'/'corr'
#ref to %iloc_cat {ilocid}='ex'/'err'/'corr'
sub perc_mis_out{
	#Declare and initialize, _r means a reference
	my($outfn,$igt_cat_r,$iloc_cat_r) = @_;
	my $niloc = 0;#total number of loci in inferred dataset
	my $nind = 0;#number of individuals
	my $ind = '';#ind id
	my $percmis = 0;
	my %ind_perc_mis = ();# {ind} percent of missing genotypes
	my $min = 0;
	my $max = 0;
	my $total = 0;
	my $avg = 0;
	my $success = 'y';
	
	unless(open(OUTFILE, ">", $outfn)){
		$success = "Cannot open outfile $outfn.";
		return $success;
	}	
	$niloc = keys %{$iloc_cat_r};
	for $ind (keys %{$igt_cat_r}){
		$percmis = $$igt_cat_r{$ind}{'mis'} / $niloc * 100;
		$ind_perc_mis{$ind} = $percmis;
	}
	$min = (sort {$a <=> $b} values %ind_perc_mis)[0];
	$max = (reverse sort {$a <=> $b} values %ind_perc_mis)[0];
	$nind = keys %{$igt_cat_r};
	for $ind (keys %ind_perc_mis){
		$total += $ind_perc_mis{$ind};
	}
	$avg = $total / $nind;
	print OUTFILE "ind\tpercent_missing_gt\n",
	"avg\t$avg\nmin\t$min\nmax\t$max\n";
	for $ind (sort keys %ind_perc_mis){
		print OUTFILE "$ind\t$ind_perc_mis{$ind}\n";
	}
	close OUTFILE;	
	return $success;
}

#Definition of subroutine perc_incorr_out
#Determines percentage of incorrect genotypes for each ind
#based on total number of loci observed in that ind
#Determines avg, min, max across inds
#Produces outfile with these values
#returns 'y' or error message as scalar
#expects
#$outfn outfilename
#ref to %igt_cat {ind}{cat}=count cat: 'mis'/'incorr'/'corr'
sub perc_incorr_out{
	#Declare and initialize, _r means a reference
	my($outfn,$igt_cat_r) = @_;
	my $ind = '';#ind id
	my $nind = 0;#number of individuals
	my $nindloc = 0;#total number of loci observed in this ind
	my $percincorr = 0;
	my %ind_perc_incorr = ();# {ind} percent of incorrect genotypes
	my $min = 0;
	my $max = 0;
	my $total = 0;
	my $avg = 0;
	my $success = 'y';
	
	unless(open(OUTFILE, ">", $outfn)){
		$success = "Cannot open outfile $outfn.";
		return $success;
	}	
	for $ind (keys %{$igt_cat_r}){
		$nindloc = $$igt_cat_r{$ind}{'incorr'} + $$igt_cat_r{$ind}{'corr'};
		if($nindloc > 0){
			$percincorr = $$igt_cat_r{$ind}{'incorr'} / $nindloc * 100;
		} else {
			$percincorr = 0;
		}
		$ind_perc_incorr{$ind} = $percincorr;
		
	}
	$min = (sort {$a <=> $b} values %ind_perc_incorr)[0];
	$max = (reverse sort {$a <=> $b} values %ind_perc_incorr)[0];
	$nind = keys %{$igt_cat_r};
	for $ind (keys %ind_perc_incorr){
		$total += $ind_perc_incorr{$ind};
	}
	$avg = $total / $nind;
	print OUTFILE "ind\tpercent_incorrect_gt\n",
	"avg\t$avg\nmin\t$min\nmax\t$max\n";
	for $ind (sort keys %ind_perc_incorr){
		print OUTFILE "$ind\t$ind_perc_incorr{$ind}\n";
	}
	close OUTFILE;	
	return $success;
}

#Definition of subroutine nloc_nind_out
#Produces outfile with number of loci with at least n individuals
#Expects:
#outfn outfilename
#ref to %iloc_nind {ilocid}=nind for each iloc: number of individuals with
				#genotypes that are finally not categorized as missing (correct and incorrect genotypes)
#returns 'y' or error message as scalar
sub nloc_nind_out{
	#Declare and initialize, _r means a reference
	my($outfn,$iloc_nind_r) = @_;
	my $ilocid = 0;
	my $nind = 0;
	my $nloc = 0;
	my %nind_nloc = ();# {nind}=nloc
	my $nindmax = 0;# greatest number of inds for any iloc
	my %atleastnind_nloc = ();# {at least n ind} = n iloc
	my $success = 'y';
	my $i = 0;
	my $j = 0;
	
	#open outfile
	unless(open(OUTFILE, ">", $outfn)){
		print "Cannot open $outfn.\n";
		$success = "Cannot open $outfn.";
		return $success;
	}
	#populate %nind_nloc
	for $ilocid (keys %{$iloc_nind_r}){
		$nind = $$iloc_nind_r{$ilocid};
		++$nind_nloc{$nind};
	}
	#get max nind for any iloc
	$nindmax = (reverse sort {$a <=> $b} keys %nind_nloc)[0];
	#initialize counts in %atleastnind_nloc with 0
	for($nind = 1; $nind <= $nindmax; ++$nind){
		$atleastnind_nloc{$nind} = 0;
	}
	#populate %atleastnind_nloc with counts
	$nloc = 0;
	for($i = 1; $i <= $nindmax; ++$i){
		for($j = $i; $j <= $nindmax; ++$j){
			if(defined $nind_nloc{$j}){
				$nloc += $nind_nloc{$j};
			}
		}
		$atleastnind_nloc{$i} = $nloc;
		$nloc = 0;
	}
	#print data to outfile
	print OUTFILE "at_least_n_ind\tn_loc\n";
	for $nind (sort {$a <=> $b} keys %atleastnind_nloc){
		print OUTFILE "$nind\t$atleastnind_nloc{$nind}\n";
	}
	close OUTFILE;
	return $success;
}
