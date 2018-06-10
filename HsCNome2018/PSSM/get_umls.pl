#!/usr/bin/perl
use strict;

print "Enter a file w UMLS IDs\n";my $idpath=<STDIN>;chomp($idpath);
open(fr,$idpath)||print "cannot open UMLS ID list $!\n";
my @data=<fr>;chomp(@data);
print "Enter the DiseaseID-DiseaseName mapping output\n";my $out=<STDIN>;chomp($out);
open(fw,">>$out");
open(fr2,"/DisGeNet_Db_v5/disease_mappings.tsv")||print "cannot open UMLS-disease_name mapping $!\n";
my @mapping=<fr2>;chomp(@mapping);my $len=@mapping;print "$len\n";#exit();

my @splitted_line;my $id;my $nomatch=0;

foreach my $line(@data)
{
 $line=~s/\n//g;$line=~s/\r//g;
 foreach my $line2(@mapping)
 {
  @splitted_line=split("\t",$line2);
  if($splitted_line[0] eq $line)
  {
   print fw "$line\t$splitted_line[1]\n";
   last();
  }
  else
  {
   $nomatch++;#print "$nomatch\t";
   if($nomatch eq $len)
   {
    print fw "$line\tNoName\n";$nomatch=0;#exit();
   }
  }
 }
 $nomatch=0;#exit();
}
