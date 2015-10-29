################################################################################################
#This scripts get all NON-overlapping putative recombinant events 
#to qualify for recombinant events, each event's recombinant gene have to match
#the RDP output.csv filename, which meant to search recombinant in that gene
#INPUT:9_seq_input_for_RDP_aligned folder (post RDP4 run, with csv files)
#OUTPUT >putative_recomb_events_NONoverlap.txt
###############################################################################################
#! /usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage; 
my $infasta='';
my $message_text  = "please specify option parameters. -i <input fasta file which contains all sequence of the family being analyzed>\n";
  my $exit_status   = 2;          ## The exit status to use
  my $verbose_level = 0;          ## The verbose level to use
  my $filehandle    = \*STDERR;   ## The filehandle to write to
GetOptions ('i=s' =>\$infasta);
pod2usage( { -message => $message_text ,
               -exitval => $exit_status  ,  
               -verbose => $verbose_level,  
               -output  => $filehandle } ) if ($infasta eq '');
			   

###process detected recombination
###get detected start and end for each event#
opendir DIR,"9_seq_input_for_RDP_aligned" or die;
open TEMPEVENT, ">putative_recomb_events_NONoverlap.txt" or die;
open TEMP,">delete_me.txt"; # store temp;
print TEMP "Recombinant gene\tMinor donor\tMajor donor\tBreakpoint start in alignment\tStart without gap\tEnd without gap\tRecombinant region length\tRDP\tGENECONV\tBootscan\tMaxchi\tChimaera\tSiSscan\tPhylPro\tLARD\tthreeSeq\n";
my %detectedrecombdetail;
my %detectedrecombdetMinor;
my %detectedrecombdetMajor;
my %recomb_num_each_gene;
my %recombMajor;#record the recomb-minor and recomb-donor combinations and start point
my %recombMinor;
my $total_recomb_event_num=0;
ONEFILE:while(my $filename=readdir(DIR))
{  
next if !($filename=~/\.csv$/);  #only scan for .csv file
next if ($filename eq "." or $filename eq ".."); 
my $gene_flag=0;
#print TEMPEVENT "\n".$filename."\n";
open DETECTEDRECOMBDETIAL, "9_seq_input_for_RDP_aligned/$filename" or die "9_seq_input_for_RDP_aligned/$filename";
my $filenamewithNo=$filename;
$filenamewithNo=~s/\.fasta\.out\.parsed\.aln\.fas\.csv//;
$filename=~s/\.fasta\.out\.parsed\.aln\.fas\.csv//;
#print $filename."\n";
my $line1=<DETECTEDRECOMBDETIAL>;
if ($line1=~/no recombination detected/ ) {  #####no recombination detected by RDP
#print "NO recombintated detected by RDP! ${filequeryname} \n";
next ONEFILE;
}
<DETECTEDRECOMBDETIAL>;
<DETECTEDRECOMBDETIAL>;#read header
my $maxAlgor=0;
my $breakflag=0;
RECOMBEVENT: while(my $line=<DETECTEDRECOMBDETIAL>) #this the the loop that scans each recomb event in a file
{
	my %event_in_current_gene; # record positions of event in current gene, use to eliminate superimposing events
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
    $RecombSeqs=~s/Unknown \(//;
	$MinorParenSeqs=~s/Unknown \(//;
	$MajParenSeqs=~s/Unknown \(//;
    $RecombSeqs=~s/REPEAT//;
	$MinorParenSeqs=~s/REPEAT//;
	$MajParenSeqs=~s/REPEAT//;
	$RecombSeqs=~s/\)//;
	$MinorParenSeqs=~s/\)//;
	$MajParenSeqs=~s/\)//;
    $RecombEnd=~s/\*//;
	$RecombBegin=~s/\*//;
	 my $recombBeginWOgap='';
	 my $recombEndWOgap='';
	#print "$RecombSeqs\t$MinorParenSeqs\t$MajParenSeqs\n";
	#print "$Maxchi\n";
	#print "detected recombseq: $RecombSeqs\n";
	my $tempRecomb='';
	my $tempMinor='';
	my $tempMaj='';
	#print "$RecombSeqs\t$filename\n";
		  if ($RecombSeqs eq $filename){
			  	
			  $AlnBegin=~s/\*//;
			  $AlnEnd=~s/\*//;
			  $AlnBegin=~s/\s//;
			  $AlnEnd=~s/\s//;
			  $recombBeginWOgap=findPositionWOgap($AlnBegin,$RecombSeqs,$filenamewithNo);
			  $recombEndWOgap=findPositionWOgap($AlnEnd,$RecombSeqs,$filenamewithNo);
			  my $recomblen=($recombEndWOgap-$recombBeginWOgap+1);
      	chomp $RecombSeqs;
	chomp $MajParenSeqs;
	chomp $MinorParenSeqs;
      if (!((exists $event_in_current_gene{'start'}{$recombBeginWOgap}) and (exists $event_in_current_gene{'end'}{$recombEndWOgap}))) {    
      print TEMP "$RecombSeqs\t$MinorParenSeqs\t$MajParenSeqs\t$AlnBegin\t$recombBeginWOgap\t$recombEndWOgap\t$recomblen\t$RDP\t$GENECONV\t$Bootscan\t$Maxchi\t$Chimaera\t$SiSscan\t$PhylPro\t$LARD\t$threeSeq\n";  ### printing current, event,the event is witin the gene we are looking for
	  $detectedrecombdetail{$RecombSeqs}{$MinorParenSeqs}{$MajParenSeqs}{$recombBeginWOgap}=1;
	  $detectedrecombdetMinor{$RecombSeqs}{$MinorParenSeqs}{$recombBeginWOgap}=($recombEndWOgap-$recombBeginWOgap+1);
	  $detectedrecombdetMajor{$RecombSeqs}{$MajParenSeqs}{$recombBeginWOgap}=($recombEndWOgap-$recombBeginWOgap+1);
        if (exists $recombMajor{$RecombSeqs}{$MajParenSeqs}{'start'}) { #record the recomb-minor and recomb-donor combinations and start point
    $recombMajor{$RecombSeqs}{$MajParenSeqs}{'start'}=$recombMajor{$RecombSeqs}{$MajParenSeqs}{'start'}."_".$recombBeginWOgap;
	}else { $recombMajor{$RecombSeqs}{$MajParenSeqs}{'start'}=$recombBeginWOgap;}
  if (exists $recombMinor{$RecombSeqs}{$MinorParenSeqs}{'start'}) { #record the recomb-minor and recomb-donor combinations
    $recombMinor{$RecombSeqs}{$MinorParenSeqs}{'start'}=$recombMinor{$RecombSeqs}{$MinorParenSeqs}{'start'}."_".$recombBeginWOgap;
	}else { $recombMinor{$RecombSeqs}{$MinorParenSeqs}{'start'}=$recombBeginWOgap;}
	  
	  if (exists $recomb_num_each_gene{$RecombSeqs}) { #record how many recomb event each gene have
	  $recomb_num_each_gene{$RecombSeqs}++;} 
	  else {$recomb_num_each_gene{$RecombSeqs}=1;}
	  $total_recomb_event_num++;
	  $event_in_current_gene{'start'}{$recombBeginWOgap}=1; #record current start and end
	  $event_in_current_gene{'end'}{$recombEndWOgap}=1;
	  }
	  }
	while (my $line2=<DETECTEDRECOMBDETIAL>)   ###search result within a recomb event(the rest of the "block in the .csv file")
		{
		next RECOMBEVENT if (!($line2=~/\d/)); # detect if it is blank line
		chomp $line2;
      my ($RecombEventNumE,$NumInRDPFileE,$AlnBeginE,$AlnEndE,$RecombBeginE,$RecombEndE,$relativeBeginE,$relativeEndE,$RecombSeqsE,$MinorParenSeqsE,$MajParenSeqsE)=split(',',$line2);

		if (($RecombSeqsE eq "") and ($tempRecomb eq '')) #no current value and no temp value, take first value of event
			{
            $RecombSeqsE=$RecombSeqs;
			}
       elsif ($RecombSeqsE eq "")     #no current value, take temp
		   {$RecombSeqsE=$tempRecomb
		   }
	     $tempRecomb=$RecombSeqsE;  #update temp

		if (($MinorParenSeqsE eq "") and ($tempMinor eq '')) #no current value and no temp value, take first value of event
			{
            $MinorParenSeqsE=$MinorParenSeqs;
			}
       elsif ($MinorParenSeqsE eq "")     #no current value, take temp
		   {$MinorParenSeqsE=$tempMinor
		   }
	     $tempMinor=$MinorParenSeqsE;  #update temp

		if (($MajParenSeqsE eq "") and ($tempMaj eq '')) #no current value and no temp value, take first value of event
			{
            $MajParenSeqsE=$MajParenSeqs;
			}
       elsif ($MajParenSeqsE eq "")     #no current value, take temp
		   {$MajParenSeqsE=$tempMaj
		   }
	     $tempMaj=$MajParenSeqsE;  #update temp  
	$RecombSeqsE=~s/Unknown\(//;
	$MinorParenSeqsE=~s/Unknown\(//;
	$MajParenSeqsE=~s/Unknown\(//;
    $RecombSeqsE=~s/REPEAT//;
	$MinorParenSeqsE=~s/REPEAT//;
	$MajParenSeqsE=~s/REPEAT//;
	$RecombSeqsE=~s/\)//;
	$MinorParenSeqsE=~s/\)//;
	$MajParenSeqsE=~s/\)//;
    $RecombEndE=~s/\*//;
	$RecombBeginE=~s/\*//;

	chomp $RecombSeqsE;
	chomp $MajParenSeqsE;
	chomp $MinorParenSeqsE;
	$MajParenSeqsE=~s/[^a-zA-Z0-9._-]+$//;
    next if (!($MajParenSeqsE =~/[a-zA-Z0-9]/) or !($MinorParenSeqsE =~/[a-zA-Z0-9]/) or !($RecombSeqsE =~/[a-zA-Z0-9]/));
		#print "$RecombSeqsE\t$MinorParenSeqsE\t$MajParenSeqsE\n";
	  if ($RecombSeqsE eq $filename){  ###the event is within the gene we are looking for
		  $AlnBegin=~s/\*//;
          $AlnEnd=~s/\*//;			
		  $AlnBegin=~s/\s//;
		  $AlnEnd=~s/\s//;
		  if ($recombBeginWOgap eq '') {
			  $recombBeginWOgap=findPositionWOgap($AlnBegin,$RecombSeqs,$filenamewithNo);
			  $recombEndWOgap=findPositionWOgap($AlnEnd,$RecombSeqs,$filenamewithNo);
			  }
		my $recomblen=$recombEndWOgap-$recombBeginWOgap+1;
 if (!((exists $event_in_current_gene{'start'}{$recombBeginWOgap}) and (exists $event_in_current_gene{'end'}{$recombEndWOgap}))) {    
 print TEMP "$RecombSeqsE\t$MinorParenSeqsE\t$MajParenSeqsE\t$AlnBegin\t$recombBeginWOgap\t$recombEndWOgap\t$recomblen\t$RDP\t$GENECONV\t$Bootscan\t$Maxchi\t$Chimaera\t$SiSscan\t$PhylPro\t$LARD\t$threeSeq\n"; 
 $detectedrecombdetail{$RecombSeqsE}{$MinorParenSeqsE}{$MajParenSeqsE}{$recombBeginWOgap}=1;
 $detectedrecombdetMinor{$RecombSeqsE}{$MinorParenSeqsE}{$recombBeginWOgap}=($recombEndWOgap-$recombBeginWOgap+1);
 $detectedrecombdetMajor{$RecombSeqsE}{$MajParenSeqsE}{$recombBeginWOgap}=($recombEndWOgap-$recombBeginWOgap+1);
        if (exists $recombMajor{$RecombSeqsE}{$MajParenSeqsE}{'start'}) { #record the recomb-minor and recomb-donor combinations and start point
    $recombMajor{$RecombSeqsE}{$MajParenSeqsE}{'start'}=$recombMajor{$RecombSeqsE}{$MajParenSeqsE}{'start'}."_".$recombBeginWOgap;
	}else { $recombMajor{$RecombSeqsE}{$MajParenSeqsE}{'start'}=$recombBeginWOgap;}
  if (exists $recombMinor{$RecombSeqsE}{$MinorParenSeqsE}{'start'}) { #record the recomb-minor and recomb-donor combinations
    $recombMinor{$RecombSeqsE}{$MinorParenSeqsE}{'start'}=$recombMinor{$RecombSeqsE}{$MinorParenSeqsE}{'start'}."_".$recombBeginWOgap;
	}else { $recombMinor{$RecombSeqsE}{$MinorParenSeqsE}{'start'}=$recombBeginWOgap;}
	  	
   	  if (exists $recomb_num_each_gene{$RecombSeqsE}) { #record how many recomb event each gene have
	  $recomb_num_each_gene{$RecombSeqsE}++;} 
	  else {$recomb_num_each_gene{$RecombSeqsE}=1;}
	  $total_recomb_event_num++;
	  $event_in_current_gene{'start'}{$recombBeginWOgap}=1; #record current start and end
	  $event_in_current_gene{'end'}{$recombEndWOgap}=1;
	  }
	  }
		}     
		}
		else {
			while (my $line2=<DETECTEDRECOMBDETIAL>)   ##skip this event because less than 2 algorithm validated 
		{
		 	next RECOMBEVENT if (!($line2=~/\d/)); # detect if it is blank line
		}
		}
}
my $RDPpositive=0;
}

close DETECTEDRECOMBDETIAL;
}
close TEMP;


###########################################
#Begin output number of events per gene
###########################################
my %eventcount;
my $total_recomb_gene_num=0;
foreach my $k1 (sort keys %recomb_num_each_gene) {
my $current_count=$recomb_num_each_gene{$k1};  #record how many event the current gene have
	if (exists $eventcount{$current_count}) {
	$eventcount{$current_count}++;
	}
	else {$eventcount{$current_count}=1;}
$total_recomb_gene_num++;
}

print TEMPEVENT "total putative recombinant gene number: $total_recomb_gene_num\ntotal putative recombinant event number:$total_recomb_event_num\n";

print TEMPEVENT "recombination event statistics:\n";
foreach my $k1 (sort {$a<=>$b}  keys %eventcount) # get and print the frequency of number of events in each gene
{
print TEMPEVENT $k1."\t".$eventcount{$k1}."\n";
}

print TEMPEVENT "GENE\trecombination event count\n";
foreach my $k1 (sort {$recomb_num_each_gene{$a} <=> $recomb_num_each_gene{$b} } keys %recomb_num_each_gene)  # print the number of events in each gene
{
print TEMPEVENT $k1."\t".$recomb_num_each_gene{$k1}."\n";
}
#############################
#Get length of all sequence
#############################
my %seq_length;
my $seqio=Bio::SeqIO->new(-format => 'Fasta', -file=>"$infasta");
while(my $seqobj=$seqio->next_seq)
{
my $id=$seqobj->id;
my $seqlength=$seqobj->length;
$seq_length{$id}=$seqlength;
}

###################################################################################
#count the occurance of annotated vs unannotated gene in recombinant/Maj/Min donor
###################################################################################
my %unannotated_count; #{recomb} or {maj} or {min} occurrance
my %annotated_count;   #{recomb} or {maj} or {min} occurrance
my %unannotated_name; # count the number of unannotated gene number, not occurrance
my %annotated_name; # count the number of annotated gene number, not occurrance
my %unannotated_name_in_recomb;
my %annotated_name_in_recomb; 
my %mixing;         #count the mixing of one seq different maj, min
$unannotated_count{'recomb'}=0;
$unannotated_count{'maj'}=0;
$unannotated_count{'min'}=0;
$annotated_count{'recomb'}=0;
$annotated_count{'maj'}=0;
$annotated_count{'min'}=0;
foreach my $k1 (sort keys %detectedrecombdetMinor) { #count recomb and min donor   $k1 is recombinant
	my $v1=$detectedrecombdetMinor{$k1};
	my %v1=%$v1;
	    if ($k1=~/unannotated/) {
		$unannotated_count{'recomb'}++;
		$unannotated_name{$k1}=1;
		$unannotated_name_in_recomb{$k1}=1;
		
			}
	else {
		$annotated_count{'recomb'}++;
		$annotated_name{$k1}=1;
        $annotated_name_in_recomb{$k1}=1;
		
			}
	foreach my $k2 (sort keys %v1) { # k1 is recomb k2 is minor
		$mixing{$k2}{'min'}=0 if !(exists $mixing{$k2}{'min'});  # initialized current recomb
		$mixing{$k2}{'minID'}='' if !(exists $mixing{$k2}{'minID'});  # initialized current recomb
	if ($k2=~/unannotated/) {
		$unannotated_count{'min'}++;
		$unannotated_name{$k2}=1;
		$mixing{$k2}{'min'}++;
		my $idNstart=$k1.'Start@'.$recombMinor{$k1}{$k2}{'start'};
		$mixing{$k2}{'minID'}.=$idNstart."__";

		}
	else {
		$annotated_count{'min'}++ ;
		$annotated_name{$k2}=1;
		$mixing{$k2}{'min'}++;
		my $idNstart=$k1.'Start@'.$recombMinor{$k1}{$k2}{'start'};
		$mixing{$k2}{'minID'}.=$idNstart."__";
	    }
	}
}
foreach my $k1 (sort keys %detectedrecombdetMajor) {
	my $v1=$detectedrecombdetMajor{$k1};
	my %v1=%$v1;
		foreach my $k2 (sort keys %v1) { # k1 is recomb k2 is major
		$mixing{$k2}{'maj'}=0 if !(exists $mixing{$k2}{'maj'});  # initialized current recomb
		$mixing{$k2}{'majID'}='' if !(exists $mixing{$k2}{'majID'});  # initialized current recomb
	if ($k2=~/unannotated/) {
		$unannotated_count{'maj'}++;
		$mixing{$k2}{'maj'}++;
		my $idNstart=$k1.'Start@'.$recombMajor{$k1}{$k2}{'start'};
		$mixing{$k2}{'majID'}.=$idNstart."__";}
	else {
		$annotated_count{'maj'}++;
		$mixing{$k2}{'maj'}++;
		my $idNstart=$k1.'Start@'.$recombMajor{$k1}{$k2}{'start'};
		$mixing{$k2}{'majID'}.=$idNstart."__";
		}
	}
}
my $tempcounter=0;
foreach my $k1 (sort keys %detectedrecombdetail) { #count recomb 
	my $v1=$detectedrecombdetail{$k1};
	my %v1=%$v1;
	foreach my $k2 (sort keys %v1) { 
		my $v2=$v1{$k2};
		my %v2=%$v2;
         foreach my $k3 (sort keys %v2) {
			 my $v3=$v2{$k3};
			 my %v3=%$v3;
			 foreach my $k4 (sort keys %v3) {	 
			 $tempcounter++;
          if (exists $mixing{$k1}{'recomb'})
			{
              $mixing{$k1}{'recomb'}++;
			}
			else {$mixing{$k1}{'recomb'}=1;}
			}
         }
	}
}
print $tempcounter;

print TEMPEVENT "frequency of annotated vs unannotated gene in recombinant/Maj/Min donor\n";
print TEMPEVENT "unannotated in:\trecombinant\tmajor donor\tminor donor:\n";
print TEMPEVENT "\t\t".$unannotated_count{'recomb'}."\t".$unannotated_count{'maj'}."\t".$unannotated_count{'min'}."\n";
print TEMPEVENT "annotated in:\trecombinant\tmajor donor\tminor donor:\n";
print TEMPEVENT "\t\t".$annotated_count{'recomb'}."\t".$annotated_count{'maj'}."\t".$annotated_count{'min'}."\n";
#count total unannotated in result and annotated in results
my $unannotated_name_count=0; # count the number of unannotated gene number, not occurrance
my $annotated_name_count=0; # count the number of annotated gene number, not occurrance
my $unannotated_name_count_recomb=0; # count the number of unannotated gene number, not occurrance
my $annotated_name_count_recomb=0; 
foreach my $k1 (sort keys %unannotated_name) {
	$unannotated_name_count++;
}
foreach my $k1 (sort keys %annotated_name) {
	$annotated_name_count++;
}
foreach my $k1 (sort keys %unannotated_name_in_recomb) {
	$unannotated_name_count_recomb++;
}
foreach my $k1 (sort keys %annotated_name_in_recomb) {
	$annotated_name_count_recomb++;
}
print TEMPEVENT "total unannotated gene in the result:$unannotated_name_count \ntotal annotated gene in the result:$annotated_name_count\n\n";
print TEMPEVENT "total unannotated gene in the recombGene:$unannotated_name_count_recomb \ntotal annotated gene in the recombGene:$annotated_name_count_recomb\n\n";
###################################################################################################################################################################
#Output details of each event, including recombination event count  Also record recomb-minor and recomb-donor combination and corresponding start for mixing output
###################################################################################################################################################################

open TEMP, "delete_me.txt";  # the details of each events is stored temporarily in this file
<TEMP>; #header
print TEMPEVENT "Recombinant gene\tMinor donor\tMajor donor\tBreakpoint start in alignment\tStart without gap\tEnd without gap\tRecombinant region length\tgene length\tRDP\tGENECONV\tBootscan\tMaxchi\tChimaera\tSiSscan\tPhylPro\tLARD\tthreeSeq\tTOTAL EVENT COUNT OF CURRENT GENE\n";
while (my $line=<TEMP>) {
	chomp $line;
    my ($RecombSeqsE,$MinorParenSeqsE,$MajParenSeqsE,$AlnBegin,$recombBeginWOgap,$recombEndWOgap,$recomblen,$RDP,$GENECONV,$Bootscan,$Maxchi,$Chimaera,$SiSscan,$PhylPro,$LARD,$threeSeq)=split("\t", $line);
    #print "$RecombSeqsE\t$MinorParenSeqsE\t$MajParenSeqsE\n";
	print TEMPEVENT "$RecombSeqsE\t$MinorParenSeqsE\t$MajParenSeqsE\t$AlnBegin\t$recombBeginWOgap\t$recombEndWOgap\t$recomblen\t".$seq_length{$RecombSeqsE}."\t$RDP\t$GENECONV\t$Bootscan\t$Maxchi\t$Chimaera\t$SiSscan\t$PhylPro\t$LARD\t$threeSeq\t$recomb_num_each_gene{$RecombSeqsE}\n";
    	
}
##########################################################################
#Output mixing (how one seq participate in different recombination events)
##########################################################################
print TEMPEVENT "\nGENE ID\tLength\tOccurances as recombinant\tOccurances as Major donor\tOccurances as Minor donor\tAs Major Donor of\t As Minor Donor of\n";
foreach my $k1 (sort keys %mixing) {
	my $majid=$mixing{$k1}{"majID"};
	my $minid=$mixing{$k1}{"minID"};
	my $recomb=$mixing{$k1}{"recomb"};
	my $maj=$mixing{$k1}{"maj"};
	my $min=$mixing{$k1}{"min"};
	$maj=0 if ((!$maj));          #deal with NULL numbers
	$min=0 if ((!$min));
	$recomb=0 if ((!$recomb));
	$majid='NA' if ((!$majid));          #deal with NULL strings
	$minid='NA' if ((!$minid));

	print TEMPEVENT $k1."\t".$seq_length{$k1}."\t".$recomb."\t".$maj."\t".$min."\t".$majid."\t".$minid."\n";
}

#####################################################################################
#defin subs
#####################################################################################
#print findPositionWOgap(200,'330_7113.t00016');  #test this function
sub findPositionWOgap # first parameter is the position with gap, second parameter is the genename, output is position without gap
{
my ($RDPalnPos,$seqname,$filenameWithNo)=@_;
my $currentseqname='';
#print '9_seq_input_for_RDP_aligned/'.$filenameWithNo.'.fasta.out.parsed.aln.fas'."\n";
if (-e '9_seq_input_for_RDP_aligned/'.$filenameWithNo.'.fasta.out.parsed.aln.fas') {
#print '9_seq_input_for_RDP_aligned/'.$filenameWithNo.'.fasta.out.parsed.aln.fas'."exists!!!!\n";
my $seqio=Bio::SeqIO->new(-format => 'Fasta', -file=>'9_seq_input_for_RDP_aligned/'.$filenameWithNo.'.fasta.out.parsed.aln.fas');
my $seqobj;
while(!($currentseqname eq $seqname) and ($seqobj=$seqio->next_seq))  #find the requested seq
{
$currentseqname=$seqobj->id;
}
if ($RDPalnPos<($seqobj->length)) {   # alignment position minus gap number equal native seq position
$_=$seqobj->subseq(1,$RDPalnPos);
my $matchcount=()=$_=~m/\-/g;
return $RDPalnPos-$matchcount;
}
else {return $RDPalnPos;}
}
else {return "error opening file";}
}

