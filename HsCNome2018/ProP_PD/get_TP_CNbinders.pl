#!/usr/bin/perl
use strict;

open(fr,"/Hs_disorderome_pepPD_Ylva/CN_ProPD_peps_12Oct")||print "cannot open CN_ProPD_peps_12Oct $!\n";
my @TP=<fr>;chomp(@TP);
open(fr2,"/Hs_disorderome_pepPD_Ylva/all_counts_12Sept")||print "cannot open all_counts_12Sept $!\n";
my @data=<fr2>;chomp(@data);

open(fw,">>/Hs_disorderome_pepPD_Ylva/all_counts_12Oct");
open(fw_tmp,">>/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_tmp");

my @splitted_line;my %peps;my $i=0;my %genes;

foreach my $line(@TP)
{
 $line=~s/\s+//g;$line=~s/\r+//g;$line=~s/\n//g;
 $peps{$line}++;
}
my $peps=keys %peps;
print "# unique peps = $peps\n";

foreach my $pep(sort keys %peps)
{
 foreach my $line(@data)
 {
  @splitted_line=split("\t",$line);
  if($pep eq $splitted_line[2])
  {
   print fw_tmp "$line\n";
  }
 }
}
close fw_tmp;

open(fr_tmp,"/Hs_disorderome_pepPD_Ylva/all_counts_12Oct_tmp")||print "Where's tmp Mister? $!\n";
my @data_tmp=<fr_tmp>;chomp(@data_tmp);

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
