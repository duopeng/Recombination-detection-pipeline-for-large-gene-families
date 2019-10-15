#! /usr/bin/perl
##################################################################################
#this script 1.splits each sequence and 2. blast against db 3.parse result, 4.get 
#hit that are >95% identity and nonself, 5.then from those select out a subset of 
#top identity hits,6 remove recent relatives, 7.and then find identity of flanking 
#sequence 
#plese specify fasta file, in which all sequence are analyzed for recombination
my $fastafile="ts_Brazil.fasta";
#specify db name, make sure there is a BLAST db in this name!!
my $dbname="ts_Brazil.fasta";
#recent diverged paralogs retained!
##################################################################################
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
system("formatdb -p F -i $dbname");
die "please make sure the director is properly set, dbname is correct, and please make sure BLAST database exists before running this script!\n" if (!(-e $dbname."\.nsq"));
# empty "1_splitted_seq" directory
if (-e "1_splitted_seq") {
	opendir DIR,"1_splitted_seq";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 1_splitted_seq/$filename");
}
closedir DIR;
}
else {mkdir("1_splitted_seq",0777) or die;}
# empty "2_splitted_seq_BLASTout" directory
if (-e "2_splitted_seq_BLASTout") {
	opendir DIR,"2_splitted_seq_BLASTout";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 2_splitted_seq_BLASTout/$filename");
}
closedir DIR;
}
else {mkdir("2_splitted_seq_BLASTout",0777) or die;}
# empty "3_splitted_seq_nonself_hit_all" directory
if (-e "3_splitted_seq_nonself_hit_all") {
	opendir DIR,"3_splitted_seq_nonself_hit_all";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 3_splitted_seq_nonself_hit_all/$filename");
}
closedir DIR;
}
else {mkdir("3_splitted_seq_nonself_hit_all",0777) or die;}
# empty "3_splitted_seq_nonself_hit_all" directory
if (-e "4_splitted_seq_nonself_hit_top") {
	opendir DIR,"4_splitted_seq_nonself_hit_top";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 4_splitted_seq_nonself_hit_top/$filename");
}
closedir DIR;
}
else {mkdir("4_splitted_seq_nonself_hit_top",0777) or die;}
# empty "5_splitted_seq_nonself_hit_top_filtered" directory
if (-e "5_splitted_seq_nonself_hit_top_filtered") {
	opendir DIR,"5_splitted_seq_nonself_hit_top_filtered";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 5_splitted_seq_nonself_hit_top_filtered/$filename");
}
closedir DIR;
}
else {mkdir("5_splitted_seq_nonself_hit_top_filtered",0777) or die;}
# empty "6_splitted_seq_nonself_hit_top_filtered_w_flankIden" directory
if (-e "6_splitted_seq_nonself_hit_top_filtered_w_flankIden") {
	opendir DIR,"6_splitted_seq_nonself_hit_top_filtered_w_flankIden";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 6_splitted_seq_nonself_hit_top_filtered_w_flankIden/$filename");
}
closedir DIR;
}
else {mkdir("6_splitted_seq_nonself_hit_top_filtered_w_flankIden",0777) or die;}
# empty "7_stat" directory
if (-e "7_stat") {
	opendir DIR,"7_stat";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 7_stat/$filename");
}
closedir DIR;
}
else {mkdir("7_stat",0777) or die;}
# empty "8_seq_input_for_RDP" directory
if (-e "8_seq_input_for_RDP") {
	opendir DIR,"8_seq_input_for_RDP";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 8_seq_input_for_RDP/$filename");
}
closedir DIR;
}
else {mkdir("8_seq_input_for_RDP",0777) or die;}
# empty "9_seq_input_for_RDP_aligned" directory
if (-e "9_seq_input_for_RDP_aligned") {
	opendir DIR,"9_seq_input_for_RDP_aligned";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("rm 9_seq_input_for_RDP_aligned/$filename");
}
closedir DIR;
}
else {mkdir("9_seq_input_for_RDP_aligned",0777) or die;}

############################################################################
#begin splitting sequence
my $windowsize=50;
my $stepsize=50;  #stepsize=$windowsize if don't want to use sliding window
die "stepsize must be less than windowsize" if $windowsize<$stepsize;
############################################################################
my $seqio=Bio::SeqIO->new(-format => 'Fasta', -file=>"${fastafile}");
while(my $seqobj=$seqio->next_seq)
{
my $id=$seqobj->id;
my $seqlength=$seqobj->length;
#print  '>'.$id."\n";
#print $seqobj->seq."\n";

# begin split one seq
open OUT, ">1_splitted_seq/${id}.fasta" or die;

my $n=0;
my $start="";
my $end="";
while($n<(int(($seqlength-$windowsize)/$stepsize)))
{
$start=$n*$stepsize+1;
$end=$n*$stepsize+$windowsize;
print OUT ">${id}"."_"."$n\n";
print OUT $seqobj->subseq($start,$end)."\n";
$n++;
}		
close OUT;
#done splitting one seq
}
#done splitting all seqs

############
#begin BLAST
############
opendir DIR,"1_splitted_seq";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system("megablast -d ${dbname} -i ./1_splitted_seq/${filename} -o 2_splitted_seq_BLASTout/${filename}.out");   
#note system("command","argu1","argu2") doest work with megablast
}
closedir DIR;
#done BLAST

#################################
#begin parsing BLASTout
#set coverage selection threshold
my $ident_threshold=95;
my $coverage_threshold=1;
#################################
opendir DIR,"2_splitted_seq_BLASTout";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";

######parsing one result#######
my %query;
open OUT, ">3_splitted_seq_nonself_hit_all/${filename}.parsed"; 
open OUT2, ">4_splitted_seq_nonself_hit_top/${filename}.parsed"; 
#print header
print OUT "query_name"."\t"."hit name"."\t"."hit->description"."\t"."percent_identity"."\t"."coverage"."\t"."start('query')"."\t"."end('query')"."\t"."start('hit')"."\t"."end('hit')"."\n";
print OUT2 "query_name"."\t"."hit name"."\t"."hit->description"."\t"."percent_identity"."\t"."coverage"."\t"."start('query')"."\t"."end('query')"."\t"."start('hit')"."\t"."end('hit')"."\n";
#打开blast结果文件
my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => "2_splitted_seq_BLASTout/${filename}");                          
my $queryno=1;  #query number, used to sort query in hash
#外层循环：扫描results
while( my $result = $in->next_result ) {
	#print OUT "next result\n";
        my $q_len = $result->query_length;
		my $queryid=$result->query_name;
		#my @queryname=split("_",$queryname);  # get the query unique ID, use later to eliminate self hit later
		#my $queryid=$queryname[1].'.'.$queryname[2];
		#print $queryid."\n";
		#print "queryid: $queryid\n";
#取all hits
  while(my $hit = $result->next_hit) {
   # print OUT "next hit\n";
    my $hit_desc=$hit->description;
    #print "hit_desc: $hit_desc\n";
#扫描hsp
   while( my $hsp = $hit->next_hsp) {   
  #  print OUT "next hsp\t";
  #  print OUT $hsp->percent_identity."\t";
   my $coverage=($hsp->length('query'))/$q_len;
   my $ident=$hsp->percent_identity;
 # print $coverage."\t";
 # print $ident."\n";
     if ($coverage>=$coverage_threshold and $ident>=$ident_threshold)     
   {
 # print OUT $result->query_name."\t".$hit->name."\t".$hit->description."\t".$hsp->start('query')."\t".$hsp->end('query')."\t".$hsp->start('hit')."\t".$hsp->end('hit')."\t".$hsp->percent_identity."\t".$coverage."\n";
    my $tempdesc=$hit->name;
	#print "hit name: $tempdesc\n";	
    if((exists $query{$queryno}) and (!($queryid=~/$tempdesc/)))
    {
     #print "query - hit  $queryid  -  $tempdesc \n";
     my $tdesc=$result->query_name."\t".$hit->name."\t".$hit->description."\t".$hsp->percent_identity."\t".$coverage."\t".$hsp->start('query')."\t".$hsp->end('query')."\t".$hsp->start('hit')."\t".$hsp->end('hit')."\n";

	 if ($ident>$query{$queryno}{'identity'}) {   #如果identity 更高，则更新
	 $query{$queryno}{'desc'}=$tdesc;	 
	 $query{$queryno}{'identity'}=$ident
	 }	 
	 if ($ident==$query{$queryno}{'identity'}) {  #如果identity 相同，则append
      $query{$queryno}{'desc'}= $query{$queryno}{'desc'}.$tdesc;	
	 }
     print OUT $tdesc;
  }
    elsif ((!(exists $query{$queryno})) and (!($queryid=~/$tempdesc/))) {  #如果没有此query的记录，则建立
          my $tdesc=$result->query_name."\t".$hit->name."\t".$hit->description."\t".$hsp->percent_identity."\t".$coverage."\t".$hsp->start('query')."\t".$hsp->end('query')."\t".$hsp->start('hit')."\t".$hsp->end('hit')."\n";
             print OUT $tdesc;
			   $query{$queryno}{'desc'}=$tdesc;
               $query{$queryno}{'identity'}=$ident;
    }
}
}
   }  
   $queryno++;
   #print $queryno."\n";
  }
  foreach my $k1 (sort {$a<=>$b} keys  %query)
{
print OUT2 $query{$k1}{'desc'};
	}
close OUT;
close OUT2;
#####end parsing one result#######
}
closedir DIR;
######done parsing


####################################################################################################################
#start to get rid of recent diverged relative sequenes and conserved region
#remove hit-seq if it max_blockno - min_blockno > 1/2 total blockno
my $total_seq_num=1446;
my $span_limit=5;#set max block span that a hit can span over the query seq
my $percent_con=0.05;#set max percent that a query block can have such many hits without being considered as conserved region
#$conserved_threshold=$percent_con*$total_seq_num;
# or set threshold for max hit to be considered as conserved region
my $conserved_threshold=30;
####################################################################################################################
#intialize some hash for storing global statistics
my %span_stat;
my $total_recomb_seq_num_unfiltered=0;
my $total_recomb_event_unfiltered=0;
my $total_recomb_seq_num_filtered=0;
my $total_recomb_event_filtered=0;
######start scanning through the folder
opendir DIR,"4_splitted_seq_nonself_hit_top";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
open INFILE, "4_splitted_seq_nonself_hit_top/$filename";
open OUTFILE, ">5_splitted_seq_nonself_hit_top_filtered/$filename";
#print OUTFILE "hit_name\t [query_block, percent_identity, left_identity, right_identity]\n";
print OUTFILE "queryname\thit_name\thit_desc\tpercent_identity\tcoverage\tstart_query\tend_query\tstart_hit\tend_hit\n";

	<INFILE>; #read header
my $blockmax=1;
my %blockstat;
my %hit_genes;
my $query_name="";
my $counter=0;
while (my $line=<INFILE>)  
{
	
		chomp $line;
		my ($querynamestring,$hit_name,$hit_desc,$percent_identity,$coverage,$start_query,$end_query,$start_hit,$end_hit)=split /\t/,$line; 
        my ($first,$VER,$second,$blockno);
		if ($querynamestring=~/VER/) { #get block No, some inconsistency in naming theme
    ($first,$VER,$second,$blockno)=split /_/, $querynamestring;
	}
	else{
    ($first,$second,$blockno)=split /_/, $querynamestring;
	 }
	  $hit_genes{$hit_name}{$blockno}{'iden'}=$percent_identity;
	  $hit_genes{$hit_name}{$blockno}{'hit_desc'}=$hit_desc;
	  $hit_genes{$hit_name}{$blockno}{'coverage'}=$coverage;
	  $hit_genes{$hit_name}{$blockno}{'start_query'}=$start_query;
	  $hit_genes{$hit_name}{$blockno}{'end_query'}=$end_query;
	  $hit_genes{$hit_name}{$blockno}{'start_hit'}=$start_hit;
	  $hit_genes{$hit_name}{$blockno}{'end_hit'}=$end_hit;
	  $hit_genes{$hit_name}{$blockno}{'querynamestring'}=$querynamestring;
	  $hit_genes{$hit_name}{$blockno}{'hit_name'}=$hit_name;
	  $blockmax=$blockno;
	  if (exists  $blockstat{$blockno}) {$blockstat{$blockno}++}
	  else {$blockstat{$blockno}=1}
	  $counter++;
     # print "$blockno  $left_iden $right_iden \n"
}
if ($counter>0) {$total_recomb_seq_num_unfiltered++}
$total_recomb_event_unfiltered=$total_recomb_event_unfiltered+$counter;
#print "$query_name ";
#print $blockstat{0}."\n";
#print $blockstat{1}."\n";
#print $blockmax."\n";
#   %hit_genes  k1=hit name k2=blockno 
###finished storing info for one file, start filtering
foreach my $k1 (sort keys %hit_genes)  #scan over each hit gene to determine if eliminates it
{
	my $v1=$hit_genes{$k1};
	my %v1=%$v1;
    my $maxkey = max keys %v1;   #max matching block no
    my $minkey = min keys %v1;   #min matching block no
#print "$k1 $minkey $maxkey \n";
my $max_span=$maxkey-$minkey+1;
if (exists $span_stat{$max_span}) {$span_stat{$max_span}++}
else{$span_stat{$max_span}=1;}

    foreach my $k2 (sort keys %v1) {
     #print $k2." ".$blockstat{$k2}."\n";
    # print $blockstat{$k2}."\t"; 
	 #print $percent_con*$total_seq_num."\n";
	 if ($blockstat{$k2}>=$conserved_threshold) {    #is conserved region if too many hits on this block
		 delete $hit_genes{$k1}{$k2};
     }
	}  
	###after going through all blocks a hit sequence have, check if all block for that hit seq is deleted
	  $v1=$hit_genes{$k1};
	 %v1=%$v1;
	  if(!%v1) {
		 # print "empty"; 
		  delete $hit_genes{$k1}};   #delete the hit sequence record if all hit blocks are deleted
 
  }


###finished filtering, starting outputting
if (!(!%hit_genes)) {           # if non-empty hit-gene hash, the sequence have recombination
	$total_recomb_seq_num_filtered++;
}
foreach my $k1 (sort keys %hit_genes)
{
	$total_recomb_event_filtered++; 
	my $v1=$hit_genes{$k1};
	my %v1=%$v1;
 foreach my $k2 (sort {$a<=>$b} keys %v1) {
	 print OUTFILE 	  $hit_genes{$k1}{$k2}{'querynamestring'}."\t".$hit_genes{$k1}{$k2}{'hit_name'}."\t".$hit_genes{$k1}{$k2}{'hit_desc'}."\t".$hit_genes{$k1}{$k2}{'iden'}."\t".$hit_genes{$k1}{$k2}{'coverage'}."\t".$hit_genes{$k1}{$k2}{'start_query'}."\t".$hit_genes{$k1}{$k2}{'end_query'}."\t".$hit_genes{$k1}{$k2}{'start_hit'}."\t".$hit_genes{$k1}{$k2}{'end_hit'}."\n";
 }
}
close INFILE;
close OUTFILE;
}

##################################################################
#start to get flanking sequence identity of hit with that of query
my $bypass_flank_iden=1;
##################################################################
my %all_seq;
############
#get all seq
############
opendir DIR,"5_splitted_seq_nonself_hit_top_filtered";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
#read in all ts seqs
my $seqio2=Bio::SeqIO->new(-format => 'Fasta', -file=>"${dbname}");
while(my $seqobj=$seqio2->next_seq)
{
my $id=$seqobj->id;
my $seq=$seqobj->seq;
$all_seq{$id}=$seq;
}
}
##########code below are for not bypass flank iden#################
if ($bypass_flank_iden==0) 
{
######start scanning through 
opendir DIR,"5_splitted_seq_nonself_hit_top_filtered";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
open INFILE, "5_splitted_seq_nonself_hit_top_filtered/$filename";
open OUTFILE, ">6_splitted_seq_nonself_hit_top_filtered_w_flankIden/$filename";
	<INFILE>; #read header
print OUTFILE "querynamestring\thit_name\thit_desc\tpercent_identity\tleft_identity\tright_identity\tcoverage\tstart_query\tend_query\tstart_hit\tend_hit\n";

while (my $line=<INFILE>)
{
		chomp $line;
		my ($querynamestring,$hit_name,$hit_desc,$percent_identity,$coverage,$start_query,$end_query,$start_hit,$end_hit)=split /\t/,$line; 
         my ($first,$VER,$second,$blockno);
		if ($querynamestring=~/VER/) { #get block No, some inconsistency in naming theme
    ($first,$VER,$second,$blockno)=split /_/, $querynamestring;
	}
	else{
    ($first,$second,$blockno)=split /_/, $querynamestring;
	 }
      my $query_name=$first."_".$second;
      #print "$query_name\n";
	  my $left_ident="NA";
	  my $right_ident="NA";
  if ((exists $all_seq{$hit_name}) and (exists $all_seq{$query_name})) {
      #print "exist";
####get query seq
####get hit len
     my $query_seq=$all_seq{$query_name};  #get full seq
     my $hit_seq=$all_seq{$hit_name};
     my $query_len=length($query_seq);     #get length
     my $hit_len=length($hit_seq);
     #print "$query_len  $hit_len  \n";
die "$query_name $hit_name" if ($start_query>$end_query or $start_hit>$end_hit);
###get left block identity
if ($blockno>0 and $start_hit>50) {     #left portion, if exists,  note block number starts at 0
my $querysubstart=($blockno-1)*50;
my $query_leftsubseq=substr($query_seq, $querysubstart,50);
my $hitsubstart=$start_hit-51;
my $hit_leftsubseq=substr($hit_seq, $hitsubstart,50);
#print "$query_leftsubseq\n$hit_leftsubseq\n";
open ALN, ">tempaln.fasta";
print ALN ">query\n$query_leftsubseq\n>hit\n$hit_leftsubseq\n";
system("clustalw2 -INFILE=tempaln.fasta -STATS=tempstat -QUIET");
open STATS, "tempstat";
for (my $i=1; $i<=13;$i++) {
	<STATS>;
}
my $line=<STATS>;
chomp $line;
my @line=split(" ",$line);
#print $line[3]."\n";  #check percent identity calculated by clustalw2
$left_ident=$line[3];
system("rm tempaln.fasta");
system("rm tempstat");
close ALN;
close STATS;
}
###get right block identity
if ((($blockno+1)*50)<($query_len-50) and $end_hit<($hit_len-50)) {     #get right portion, if exists
my $querysubstart=($blockno+1)*50;
my $query_rightsubseq=substr($query_seq, $querysubstart,50);
my $hitsubstart=$end_hit;
my $hit_rightsubseq=substr($hit_seq, $hitsubstart,50);
#print "$query_leftsubseq\n$hit_leftsubseq\n";
open ALN, ">tempaln.fasta";
print ALN ">query\n$query_rightsubseq\n>hit\n$hit_rightsubseq\n";
system("clustalw2 -INFILE=tempaln.fasta -STATS=tempstat -QUIET");
open STATS, "tempstat";
for (my $i=1; $i<=13;$i++) {
	<STATS>;
}
my $line=<STATS>;
chomp $line;
my @line=split(" ",$line);
#print $line[3]."\n";  #check percent identity calculated by clustalw2
$right_ident=$line[3];
system("rm tempaln.fasta");
system("rm tempstat");
close ALN;
close STATS;
}
	}
print OUTFILE $querynamestring."\t".$hit_name."\t".$hit_desc."\t".$percent_identity."\t".$left_ident."\t".$right_ident."\t".$coverage."\t".$start_query."\t".$end_query."\t".$start_hit."\t".$end_hit."\n";
	}
close INFILE;
close OUTFILE;
#print "done $filename \n";
}
close DIR;
#######finished getting flanking sequence identity
}

########################################################################################################
#start calculating statistics
my $histo_interval=5;
my $inter_num=60;  #number of intervals
#other needed stat: total retained hits, total seqs have hits, spanning length(get from simulated data)
########################################################################################################
if ($bypass_flank_iden==1) 
{
open OUT, ">7_stat/hit_number_count.txt";
my %blockhisto;
my %blockhistogram;
for (my $i=1;$i<=$inter_num;$i++) {   #initialized hash to store final histgram values
$blockhistogram{$i*$histo_interval}=0;
}

opendir DIR,"5_splitted_seq_nonself_hit_top_filtered";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
open INFILE, "5_splitted_seq_nonself_hit_top_filtered/$filename";
<INFILE>; #read header
my %blockstat;
while (my $line=<INFILE>)  #read through current file
{
		chomp $line;
		my ($querynamestring,$hit_name,$hit_desc,$percent_identity,$coverage,$start_query,$end_query,$start_hit,$end_hit)=split /\t/,$line; 
        my ($first,$second,$blockno)=split /_/, $querynamestring;
        if (exists  $blockstat{$blockno}) {$blockstat{$blockno}++}
	    else {$blockstat{$blockno}=1}

}
foreach my $k1 (sort keys %blockstat)  
{
  my $count=$blockstat{$k1};
  if (exists $blockhisto{$count}) {$blockhisto{$count}+=$count}
  else {$blockhisto{$count}=$count}
}
close INFILE
}
close DIR;

###start to get histogram values
 foreach my $k1 (sort {$a<=>$b} keys %blockhisto) 
{
 print OUT "blocks that matched with $k1 times have a total hit count of $blockhisto{$k1} \n";
for (my $i=1;$i<=$inter_num;$i++) 
  {
	my $uplimit=$i*$histo_interval;
	my $lowlimit=$i*$histo_interval-$histo_interval;
if (($lowlimit<$k1 and $k1<=$uplimit) or ($i==1 and $k1==1)) 	
	{
	if (exists $blockhistogram{$uplimit}) {$blockhistogram{$uplimit}+=$blockhisto{$k1}}
	else {$blockhistogram{$uplimit}=$blockhisto{$k1}}
	}
  }

}
 foreach my $k1 (sort {$a<=>$b} keys %blockhisto) 
{
 print OUT "$k1\t$blockhisto{$k1} \n";
}
###start print out histogram values
foreach my $k1 (sort {$a<=>$b} keys %blockhistogram) 
{
	my $lower=$k1-$histo_interval+1;
 print OUT "$lower-$k1\t$blockhistogram{$k1} \n";

}


###start print out max span values
print OUT "MAXIMUM SPAN\n span_length\t count\n";
foreach my $k1 (sort {$a<=>$b} keys %span_stat) 
{
 print OUT "${k1}\t".$span_stat{$k1}."\n";

}
###start print out event counts
print OUT "total_recomb_seq_num_unfiltered:\t$total_recomb_seq_num_unfiltered \n";
print OUT "total_recomb_event_unfiltered:\t$total_recomb_event_unfiltered \n";
print OUT "total_recomb_seq_num_filtered:\t$total_recomb_seq_num_filtered \n";
print OUT "total_recomb_event_filtered:\t$total_recomb_event_filtered \n";
close OUT;
}
else{  ## flank ident not bypassed

}

################################################################################################
#get potential recombinat gene groups for RDP input
################################################################################################
opendir DIR,"5_splitted_seq_nonself_hit_top_filtered";
SCANCURRENTFILE: while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
close OUTFILE;
close INFILE;
open OUTFILE, ">8_seq_input_for_RDP/$filename";
open INFILE, "5_splitted_seq_nonself_hit_top_filtered/$filename";
<INFILE>; #read header
my %blockstat;
my $queryname;
my %hitnames;
my $templine;
if (!($templine=<INFILE>)) {    # deal with empty file with only 1 line header
	system ("rm 8_seq_input_for_RDP/$filename");
	next SCANCURRENTFILE;
		}

seek(INFILE, -length($templine), 1); # place the same line back onto the filehandle
 while (my $line=<INFILE>)  #read through current file
{   
		chomp $line;
		my ($querynamestring,$hit_name,$hit_desc,$percent_identity,$coverage,$start_query,$end_query,$start_hit,$end_hit)=split /\t/,$line; 
        my ($first,$VER,$second,$blockno);
		$VER="";
		if ($querynamestring=~/VER/) { #get block No, some inconsistency in naming theme
    ($first,$VER,$second,$blockno)=split /_/, $querynamestring;
	}
	else{
    ($first,$second,$blockno)=split /_/, $querynamestring;
	 }
	 if ($VER eq "") {
		    $queryname=$first."_".$second;
	 }
	 else{
		  $queryname=$first."_".$VER."_".$second;
	 }
       
  	 $hitnames{$hit_name}=1;
}
#get query seq
$queryname=~s/_.+$//;
print OUTFILE ">${queryname}\n";
print OUTFILE $all_seq{$queryname}."\n";
#get hit seqs
foreach my $k1 (sort keys %hitnames) {
print OUTFILE ">${k1}\n";
print OUTFILE $all_seq{$k1}."\n";
}
}

#####################################################################
#start align groups of potential recombinant sequences
#####################################################################
opendir DIR,"8_seq_input_for_RDP";
while(my $filename=readdir(DIR))
{  
next if $filename eq "." or $filename eq "..";
system ("clustalw2 8_seq_input_for_RDP/$filename -ALIGN -OUTPUT=fasta -OUTORDER=INPUT -QUIET -OUTFILE=9_seq_input_for_RDP_aligned/$filename.aln.fas")
}
