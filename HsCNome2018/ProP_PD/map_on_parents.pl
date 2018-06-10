#!/usr/bin/perl
use strict;

#open(fr,"/Hs_disorderome_pepPD_Ylva/all_counts_12Oct")||print "cannot open all_counts_12Oct $!\n";          # Round 1
open(fr,"/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_corrected")||print "cannot open all_counts_12Oct $!\n"; # Round 2
my @data=<fr>;chomp(@data);
open(fw_tmp,">>/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_tmp");
open(fw2,">>/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_peps_unmapped_round2");
my @splitted_line;my %UPs;my $end;my $start;my @data2;my %genes;

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);#print "$splitted_line[1]\t$splitted_line[2]\n";exit();
 $UPs{$splitted_line[1]}{$splitted_line[2]}="";
}
my $UPs=keys %UPs;print "I will search parents of $UPs Acc numbers\n";

my $dir="/Hs_proteome_10July_2K15_rev_n_unreviewed";
opendir(DIR,"$dir")||print "Cannot open parents' house $!\n";

foreach my $UP(sort keys %UPs)
{
 open(fr,"$dir/$UP\.fasta")||print "Cannot access the parent $!\n";
 @data2=<fr>;chomp(@data2);

 foreach my $pep(sort keys %{$UPs{$UP}})
 {
  if($data2[1]=~/$pep/g)
  {
   $end=pos($data2[1]);
   $start=$end - length($pep)+1;
   print "$pep matches $data2[0] i.e. $UP from $start to $end\n";

   foreach my $line(@data)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[1] eq $UP)&&($splitted_line[2] eq $pep))
    {
     print fw_tmp "$splitted_line[0]\t$splitted_line[1]\t$splitted_line[2]\t$start\t$end\t$splitted_line[3]\t$splitted_line[4]\t$splitted_line[5]\t$splitted_line[6]\t$splitted_line[7]\t$splitted_line[8]\n";
    }
   }
  }
  else
  {
   #print "I could not map $pep to $data2[0]\n";exit();
   print fw2 "$UP\t$pep\n";
  }
 }
}
close fw_tmp;

open(fr_tmp,"/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_tmp")||print "Where's tmp Mister? $!\n";
my @data_tmp=<fr_tmp>;chomp(@data_tmp);
open(fw,">>/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_with_pos_round2");
foreach my $line(@data_tmp)
{
 @splitted_line=split("\t",$line);
 $genes{$splitted_line[0]}="";
}
my $genes=keys %genes;
print "There are $genes unique genes\n";

foreach my $gene(sort keys %genes)
{
 foreach my $line(@data_tmp)
 {
  @splitted_line=split("\t",$line);
  if($gene eq $splitted_line[0])
  {
   print fw "$line\n";
  }
 }
}
system("rm -f /Hs_disorderome_pepPD_Ylva/all_counts_12Oct_tmp");
