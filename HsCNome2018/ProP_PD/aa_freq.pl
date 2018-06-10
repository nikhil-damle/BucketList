#### This prog computes amino acid freq in the disordered regions in human proteome. ####

#!/usr/bin/perl
use strict;

open(fr,"/hDisorderome/hDisorderome.txt")||print "cannot open hDisorderome.txt $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;my $freq;my $percent_freq;my %aa;my @splitted_pep;my $len;my $total_aa;

foreach my $line(@data)
{
 if($line!~/^acc\t/)
 {
  @splitted_line=split("\t",$line);
  if($splitted_line[1]!~/[BJOUXZ]/g)
  {
   $splitted_line[1]=~s/\s+//g;$splitted_line[1]=~s/\n//g;$splitted_line[1]=~s/\r//g;$splitted_line[1]=~s/[^A-Z]//g;
   @splitted_pep=split("",$splitted_line[1]);
   foreach my $aa(@splitted_pep)
   {
    $aa{$aa}++;
   }
   $len=$len+scalar @splitted_pep;
  }
 }
}
foreach my $aa(sort keys %aa)
{
 $total_aa = $total_aa+$aa{$aa};
}
print "Len = $len\nTotal aas = $total_aa\n";

foreach my $aa(sort {$aa{$b}<=>$aa{$a}} keys %aa)
{
 $freq=($aa{$aa}/$total_aa)*100;
 $freq=sprintf("%.5f",$freq);
 print "$aa\t$aa{$aa}\t$total_aa\t$freq\n";
}
