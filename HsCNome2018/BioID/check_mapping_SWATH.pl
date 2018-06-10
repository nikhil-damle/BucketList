#!/usr/bin/perl
use strict;

my $dir="SWATH";
open(fr,"$dir/count_matrix_anyBFDR_le0.05")||print "cannot open mapping file $!\n";
my @data=<fr>;chomp(@data);
open(fr2,"$dir/geneID_UP_mapping_BioMart_Ensembl89"||print "cannot open mapping $!\n");   # This file was obtained from Ensembl/BioMart at the time
my @mapping=<fr2>;chomp(@mapping);

open(fw,">>$dir/count_matrix_anyBFDR_le0.05_w_mapping");

my @splitted_line;my %genes;my %maps;my $nomatch=0;my $k;

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);
 $genes{$splitted_line[0]}="";
}
my $genes=keys %genes;print "$genes appeared w BFDR<=0.05 in at least one of the Cn expt\n";

foreach my $line(@mapping)
{
 @splitted_line=split("\t",$line);
 $maps{$splitted_line[0]}{$splitted_line[1]}="";
}
my $mapped=keys %maps;print "$mapped genes mapped onto UP acc\n";

foreach my $gene1(sort keys %genes)
{
 foreach my $line(@data)
 {
  @splitted_line=split("\t",$line);
  if($gene1 eq $splitted_line[0])
  {
   print fw "$gene1\t";
   foreach my $gene2(sort keys %maps)
   {
    if($gene1 eq $gene2)
    {
     foreach my $acc(sort keys %{$maps{$gene1}})
     {
      $acc=~s/\s+//g;$acc=~s/\t//g;
      if($acc ne ""){print fw "$acc, ";}
     }
     print fw "\t";
    }
    else
    {
     $nomatch++;
     if($nomatch eq $mapped)
     {
      print "Not mapped\t";#$gene1 could not be mapped onto UP acc, check it\n";
     }
    }
   }
   for($k=1;$k<scalar @splitted_line;$k++)
   {
    print fw "$splitted_line[$k]\t";
   }
   print fw "\n";
  }
 }
 $nomatch=0;
}
close fw;

open(fr,"$dir/count_matrix_anyBFDR_le0.05_w_mapping")||print "cannot open mapping file $!\n";
my @data=<fr>;chomp(@data);
my $i;my %hash;my %unmapped;my @splitted_line;my $count=0;
for($i=3;$i<scalar @data;$i++)
{
 @splitted_line=split("\t",$data[$i]);
 $hash{$splitted_line[0]}++;
 if($splitted_line[1] eq "Not mapped")
 {
  $unmapped{$splitted_line[0]}++;
 }
}
my $hash=keys %hash;print "You have $hash prots\n";
my $unmapped=keys %unmapped;print "$unmapped of which could not be mapped\n";

foreach my $acc(keys %hash)
{
 print "$acc\t$hash{$acc}\n";
 $count++;
 if($count eq 10)
 {
  last();
 }
}
