#! /usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SearchIO;
use List::Util 'max';
use List::Util 'min';
#############################################################################################################################
#start run RDP
#before running,DO the following:
#copy 6_splitted_seq_nonself_hit_top_filtered_w_flankIden to working DIR in windows
#copy 5_splitted_seq_nonself_hit_top_filtered to working DIR in windows
#copy 9_seq_input_for_RDP_aligned to working DIR in windows
#copy 7_stat DIR in windows
#set input directory  NOTE: no space in URL!!
my $workingdir="C:\\bioinfo\\mucin_win50";
my $infiledir=$workingdir."\\9_seq_input_for_RDP_aligned\\"; 
#set RDP directory
my $rdpURL="C:\\RDP4\\RDP4.exe";  
#OUTPUT:
#DETECTEDRECOMB.txt: significant recombination record
#RDPfail.record.txt: files that failed to analyze by RDP
#${workingdir}\\7_stat\\RDP.txt": statistics
##############################################################################################################################
my $insuff_seq_count=0;  # records the number of files that do not have enough sequence for RDP4 analysis(3 minimum)
my $no_RDP_recomb=0;
my $total_count=0;
my $RDPpositive=0;
my $geneRecombPositive=0;
my $RDPfail=0;
my %RDPresult;

open DETECTEDRECOMB, ">DETECTEDRECOMB.txt";
open RDPFAIL, ">RDPfail.record.txt";
opendir DIR,"9_seq_input_for_RDP_aligned" or die;
system ("del $infiledir\\\*\.rdp"); #remove previous generated rdp file
system ("del $infiledir\\\*\.csv"); #remove previous generated csv file
#system ("9_seq_input_for_RDP_aligned\\del \*\.csv");
ONEFILE: while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
#start count seq number
my $seqio=Bio::SeqIO->new(-format => 'Fasta', -file=>"9_seq_input_for_RDP_aligned\\$filename");
my $seqcount=0;
$total_count++;
while(my $seqobj=$seqio->next_seq)   #### check how many seq is in the file
	{
$seqcount++;
	}
	if ($seqcount>2) {   ####only run RDP if more than 2 seqs is in the file####
		print "\n${total_count}. excuting RPD on $filename\n";
		system ("$rdpURL -f $infiledir$filename -o");
		print "system return status $?  $!\n";
       if ($?==0 and (-e "9_seq_input_for_RDP_aligned\\${filename}.csv")) { #### check if RDP exucted successfully and output CSV file exists
    #### excute RDP4
	
	    #######start analysis of ONE RDP4 result ######
open RDPRESULT, "9_seq_input_for_RDP_aligned\\${filename}.csv";
my @file=split(/\./,$filename);
my $filequeryname=$file[0].'.'.$file[1];
#print $filequeryname;
my $line1=<RDPRESULT>;
if ($line1=~/no recombination detected/ ) {  #####no recombination detected by RDP
print "NO recombintated detected by RDP! ${filequeryname} \n";
$no_RDP_recomb++;
next ONEFILE;
}
$RDPpositive++;
<RDPRESULT>;
<RDPRESULT>;#read header
my $retainFlag=0;
my $maxAlgor=0;
my $breakflag=0;
RECOMBEVENT: while(my $line=<RDPRESULT>) #this the the loop that scans each recomb event in a file
{
	chomp $line;
    if ($line=~/The actual breakpoint position is undetermined/) { #check if it is end of file 
	$breakflag=1;
    last RECOMBEVENT;
	}
	if ($breakflag==0) {
	
	next RECOMBEVENT if (!($line=~/\d/)); # detect if it is blank line
	my ($RecombEventNum,$NumInRDPFile,$AlnBegin,$AlnEnd,$RecombBegin,$RecombEnd,$relativeBegin,$relativeEnd,$RecombSeqs,$MinorParenSeqs,$MajParenSeqs,$RDP,$GENECONV,$Bootscan,$Maxchi,$Chimaera,$SiSscan,$PhylPro,$LARD,$threeSeq)=split(',',$line);     
#print $line."\n";
##retain this result if >2 method detect recombination
my $catresult=$RDP.$GENECONV.$Bootscan.$Maxchi.$Chimaera.$SiSscan.$PhylPro.$LARD.$threeSeq;
#print "$catresult\n";
my $sigCount=0;
$sigCount++ while ($catresult=~m/E-/g);  # count the number of method return positive result
$maxAlgor=$sigCount if $sigCount>$maxAlgor;
if ($sigCount>2) {  #####current recomb event have >2 algorithm validation, go through all subevent to see if the query gene is recombinant
	#print "$catresult match more than 2 times\n";
    
	#print "$Maxchi\n";
	while (my $line2=<RDPRESULT>)   ##result within a recomb event
		{
		next RECOMBEVENT if (!($line2=~/\d/)); # detect if it is blank line
		chomp $line2;
      ($RecombEventNum,$NumInRDPFile,$AlnBegin,$AlnEnd,$RecombBegin,$RecombEnd,$relativeBegin,$relativeEnd,$RecombSeqs,$MinorParenSeqs,$MajParenSeqs)=split(',',$line2);
		$retainFlag++ if ($RecombSeqs=~/$filequeryname/ or $MinorParenSeqs=~/$filequeryname/ or $MajParenSeqs=~/$filequeryname/);
		#print $line2."\n"
		}

		}
		else {
			while (my $line2=<RDPRESULT>)   ##result within a recomb event
		{
		 	next RECOMBEVENT if (!($line2=~/\d/)); # detect if it is blank line
		}
		}
}
my $RDPpositive=0;
}
$geneRecombPositive++ if $retainFlag>=1;
print DETECTEDRECOMB "$retainFlag recombinations, at least one validated by $maxAlgor algorithms, $filequeryname\n" if $retainFlag>=1;
print "$retainFlag recombinations, at least one validated by $maxAlgor algorithms, $filequeryname\n";
$RDPresult{$filequeryname}{'recombNO'}=$retainFlag;
$RDPresult{$filequeryname}{'maxAlgor'}=$maxAlgor;
close RDPRESULT;
                                ######END analysis of one RDP rsult######
	}
	else {
		$RDPfail++;
		print RDPFAIL "$filename is unable to process through RDP\n";
	}
#####finsish analyzing one result
	}
	else{$insuff_seq_count++;}
}
closedir DIR;
#######start output result
open OUTRDP,">${workingdir}\\7_stat\\RDP.txt";
print OUTRDP "gene_id\trecombination_event_number\tnumber_of_algorithms_validated\n";
my $recobGeneCount=0;
my $recombEveCount=0;
foreach my $key (sort keys %RDPresult) {
	print OUTRDP $key."\t".$RDPresult{$key}{'recombNO'}."\t".$RDPresult{$key}{'maxAlgor'}."\n";
 if ($RDPresult{$key}{'recombNO'}>=1 and $RDPresult{$key}{'maxAlgor'}>=2) {
 $recobGeneCount++;
 $recombEveCount+=$RDPresult{$key}{'recombNO'};
 }
}
print OUTRDP "validated gene number: $recobGeneCount\n";
print OUTRDP "validated recombNO: $recombEveCount\n";
print OUTRDP "count that RDP detect no recomb: ${no_RDP_recomb}\n";
print OUTRDP "number of files that have less than 2 seqs(didn't do RDP analysis): ${insuff_seq_count}\n";
print OUTRDP "number of files that crashed RDP analyzes: ${RDPfail}\n";
print OUTRDP "number of files that showed positive result in RDP ${RDPpositive}\n";
print OUTRDP "number of files that showed true positive RDP result ${geneRecombPositive}\n";
print OUTRDP "number of total input files: ${total_count}";
close OUTRDP;