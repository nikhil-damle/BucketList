#!/usr/bin/perl
use strict;

open(fr,"HsCNome_UPacc_NCBIGeneID_mapping")||print "cannot open HsCNome_UPacc_NCBIGeneID_mapping_1Mar2018 $!\n";
my @data=<fr>;chomp(@data);

open(fr2,"hiConf_HsCNome_acc")||print "cannot open original acclist $!\n";   # List of accession IDs of HsCNome
my @list=<fr2>;chomp(@list);

my @splitted_line;my $i;my %UPs;my %GIDs;my %mappings;my $pair;my $nomatch=0;my %hash;

foreach my $line(@list)
{
 $line=~s/\n//g;$line=~s/\r//g;
 $hash{$line}="";
}
my $cnome=keys %hash;print "You ori list had $cnome acc\n";

for($i=1;$i<scalar @data;$i++)
{
 @splitted_line=split("\t",$data[$i]);
 if($splitted_line[0] ne ""){$UPs{$splitted_line[0]}++;}
 if($splitted_line[1] ne ""){$GIDs{$splitted_line[1]}++;}
 $pair=$splitted_line[0]."_".$splitted_line[1];
 $mappings{$pair}++;

 if(($splitted_line[0] eq "")||($splitted_line[1] eq ""))
 {
  print "check this pair: $pair\n";
 }
}
my $UPs=keys %UPs;print "You have mapped $UPs UPs ";
my $GIDs=keys %GIDs;print "to $GIDs geneIDs\n";
print "Following acc did not get mapped\n";
foreach my $acc1(sort keys %hash)
{
 foreach my $acc2(sort keys %UPs)
 {
  if($acc1 ne $acc2)
  {
   $nomatch++;
   if($nomatch eq $UPs)
   {
    print "$acc1\n";
   }
  }
 }
 $nomatch=0;
}
print "Following UP acc were mapped to multiple geneIDs\n";
foreach my $acc(sort keys %UPs)
{
 if($UPs{$acc}>1){print "$acc\n";}
}
print "Following geneIDs were mapped to multiple UPacc\n";
foreach my $id(sort keys %GIDs)
{
 if($GIDs{$id}>1){print "$id\n";}
}
