#!/usr/bin/perl
use strict;

open(fr,"top5K_hits")||print "cannot open top5K_hits $!\n";
my @data=<fr>;chomp(@data);
open(fw,">>top5K_hits_filtered");

my @splitted_line;my %motifs;my $str;

foreach my $line(@data)
{
 if($line!~/^pssm_p/)
 {
  @splitted_line=split("\t",$line);
  $str=$splitted_line[5]."_".$splitted_line[8]."_".$splitted_line[9]."_".$splitted_line[10]."_".$splitted_line[12];
  $motifs{$str}="";
 }
}
my $motifs=keys %motifs;
print "PSSM predicted $motifs motifs\n";

foreach my $line(@data)
{
 if($line!~/^pssm_p/)
 {
  @splitted_line=split("\t",$line);
  $str=$splitted_line[5]."_".$splitted_line[8]."_".$splitted_line[9]."_".$splitted_line[10]."_".$splitted_line[12];
  print fw "$line\t";
  foreach my $motif(sort keys %motifs)
  {
   if($motif eq $str)
   {
    if(($splitted_line[3]!~/Secreted/)&&($splitted_line[3]!~/Extracellular/)&&($splitted_line[3]!~/Mitochondrion/)&&($splitted_line[3]!~/Vesicular/)&&($splitted_line[3]!~/Lumenal/)&&($splitted_line[3]!~/Mitochondrial intermembrane/)&&($splitted_line[3]!~/Lysosome/)&&($splitted_line[3]!~/Lysosome/)&&($splitted_line[3]!~/Endosome/)&&($splitted_line[3]!~/transmembrane region/)&&($splitted_line[3]!~/Peroxisome/)&&($splitted_line[3]!~/Perinuclear space/)&&($splitted_line[3]!~/intramembrane region/)&&($splitted_line[3]!~/Intragranular/)&&($splitted_line[3]!~/Vacuolar/)&&($splitted_line[3]!~/Mitochondrial matrix/)&&($splitted_line[3]!~/Exoplasmic loop/))
    {
     print fw "Y\t";
    }
    else{print fw "N\t";}
    if($splitted_line[38] >= 0.3)
    {
     print fw "Y\t";
    }
    else{print fw "N\t";}
    if($splitted_line[0] <= 0.0001)
    {
     print fw "Y\n";
    }
    else{print fw "N\n";}
   }
  }
 }#exit();
}
