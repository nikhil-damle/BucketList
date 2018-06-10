#!/usr/bin/perl
use strict;

print "Enter the entire name incl of file path where you wish to search for comboSites\n";my $file=<STDIN>;chomp($file);
open(fr,$file)||print "cannot open $file $!\n";
my @data=<fr>;chomp(@data);

system("rm -f combo_sites");
open(fw,">>combo_sites");

my @splitted_line;my %prots;my $i;my $str;my $j;my $combo;

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);
 $prots{$splitted_line[5]}="";
 if($splitted_line[10]=~/L[A-Z]VP[A-Z][ILVF][A-Z][ILVF][TSHDEQNKR]/)
 {
  for($i=0;$i<12;$i++){print fw "$splitted_line[$i]\t";}
  print fw "$splitted_line[12], LxVPxIxIT\t";
  for($j=13;$j<scalar @splitted_line;$j++){print fw "$splitted_line[$j]\t";}print fw "\n";
 }
 if($splitted_line[12] eq "LxVP")
 {
  #print "Predicted $splitted_line[8]-$splitted_line[10]-$splitted_line[9]\n";
  $str=$splitted_line[8]."_".$splitted_line[9];
  #print "Passing $splitted_line[5] and $str\n";
  $combo=&get_stretch($splitted_line[5],$str);
  if($combo ne "")
  {
   for($i=0;$i<10;$i++){print fw "$splitted_line[$i]\t";}
   print fw "$splitted_line[10], $combo\t$splitted_line[11]\t$splitted_line[12], LxVPxIxIT\t";
   for($j=13;$j<scalar @splitted_line;$j++){print fw "$splitted_line[$j]\t";}print fw "\n";#exit();
  }
 }
}

sub get_stretch()
{
 my $acc=$_[0];my $str=$_[1];
 $str=~/(\d+)\_(\d+)/;
 my $start=$1;my $end=$2+4;my $end2=$2;           # $end=$2+4; end includes LxVPx part, so u need is IxIT part
 my $len=$end-$start+1;my $len2=$end2-$start+1;
 #print "$acc\t$start\t$end\t$len\n";exit()
 
 open(fr,"/features_combined/$acc\.fasta")||print "cannot open $acc\.fasta $!\n";
 my @data=<fr>;chomp(@data);

 my $motif_ori=substr($data[1],$start-1,$len2);
 my $motif=substr($data[1],$start-1,$len);
 #print "Original Predicted Motif=$start-$motif_ori-$end2\nNew Motif=$start-$motif-$end\n";exit();
 if($motif=~/L[A-Z]VP[A-Z][ILVF][A-Z][ILVF][TSHDEQNKR]/)
 {#print "$motif\n";
  return $motif;
 }
}
