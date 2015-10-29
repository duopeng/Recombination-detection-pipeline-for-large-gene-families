#####################################################################################################################
#This scripts realigns the final recombination events to remove potential alignment error while aligning the groups
my $fastafile="";
#####################################################################################################################
#! /usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage; 
my $message_text  = "please specify option parameters. -i <fasta file containing all sequences in the gene family>\n";
  my $exit_status   = 2;          ## The exit status to use
  my $verbose_level = 0;          ## The verbose level to use
  my $filehandle    = \*STDERR;   ## The filehandle to write to
GetOptions ('i=s' =>\$fastafile);
pod2usage( { -message => $message_text ,
               -exitval => $exit_status  ,  
               -verbose => $verbose_level,  
               -output  => $filehandle } ) if ($fastafile eq '');

############
#get all seq
############
my %all_seq;
my $seqio2=Bio::SeqIO->new(-format => 'Fasta', -file=>"$fastafile");
while(my $seqobj=$seqio2->next_seq)
{
my $id=$seqobj->id;
my $seq=$seqobj->seq;
$all_seq{$id}=$seq;
}

#######################################
#Get all the triplets into fas files
#######################################
open IN, "putative_recomb_events_NONoverlap.txt" or die;
while (my $line=<IN>) {  #skip till headers

	last if $line=~m/Recombinant gene	Minor donor	Major donor	Breakpoint start in alignment	Start without gap	End without gap	Recombinant region length	gene length	RDP	GENECONV	Bootscan	Maxchi	Chimaera	SiSscan	PhylPro	LARD	threeSeq	TOTAL EVENT COUNT OF CURRENT GENE/;

}
my $count=1;

if (-e "10_realign_W_triplet") {
	opendir DIR,"10_realign_W_triplet";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 10_realign_W_triplet/$filename");
}
closedir DIR;
}
else {mkdir("10_realign_W_triplet",0777) or die;}

opendir DIR,"10_realign_W_triplet";
while (my $line=<IN>) {
	#print "test\n";
	my $flag=($line=~m/GENE ID	Length	Occurances as recombinant	Occurances as Major donor	Occurances as Minor donor	As Major Donor of	 As Minor Donor of/);

	last if $flag==1;
	if ($flag!=1) {
    my ($RecombSeqsE,$MinorParenSeqsE,$MajParenSeqsE,$AlnBegin,$recombBeginWOgap,$recombEndWOgap,$recomblen,$RDP,$GENECONV,$Bootscan,$Maxchi,$Chimaera,$SiSscan,$PhylPro,$LARD,$threeSeq)=split("\t", $line);
#print "$RecombSeqsE,$MinorParenSeqsE,$MajParenSeqsE\n";
	#open the output fasta file and write the sequence
	if ($MajParenSeqsE ne '') {

###deal with same names in recomb, maj, min, remember to trim off "REPEAT"!!!!!!!!!!!;
if ($RecombSeqsE eq $MinorParenSeqsE) {$MinorParenSeqsE=$MinorParenSeqsE."REPEAT";}
if ($RecombSeqsE eq $MajParenSeqsE) {$MajParenSeqsE=$MajParenSeqsE."REPEAT";}
if ($MinorParenSeqsE eq $MajParenSeqsE) {$MajParenSeqsE=$MajParenSeqsE."REPEAT";}
	
	open OUTFILE, ">10_realign_W_triplet/${RecombSeqsE}\.No${count}\.fasta\n";
	print OUTFILE ">${RecombSeqsE}\n";
    $RecombSeqsE=~s/REPEAT//;
	print OUTFILE $all_seq{$RecombSeqsE}."\n";

	print OUTFILE ">${MajParenSeqsE}\n";
	$MajParenSeqsE=~s/REPEAT//;
	print OUTFILE $all_seq{$MajParenSeqsE}."\n";

	print OUTFILE ">${MinorParenSeqsE}\n";
	$MinorParenSeqsE=~s/REPEAT//;
	print OUTFILE $all_seq{$MinorParenSeqsE}."\n";
$count++;
}
}
}
closedir DIR;
#####################################################################
#start align groups of triplets
#####################################################################
if (-e "11_tripletALN_input_for_RDP") {
	opendir DIR,"11_tripletALN_input_for_RDP";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 11_tripletALN_input_for_RDP/$filename");
}
closedir DIR;
}
else {mkdir("11_tripletALN_input_for_RDP",0777) or die;}


opendir DIR,"10_realign_W_triplet";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
print "aligning 10_realign_W_triplet/$filename\n";
system ("clustalw2 10_realign_W_triplet/$filename -ALIGN -OUTPUT=fasta -OUTORDER=INPUT -QUIET -OUTFILE=11_tripletALN_input_for_RDP/$filename.aln.fas")
}
