#!/usr/bin/perl
use strict;

#my $dir="SWATH";
my $dir="DDA";

#open(fr,"$dir/count_matrix_SWATH_irresp_o_BFDR")||print "cannot open file $!\n";
open(fr,"$dir/count_matrix_DDA_irresp_o_BFDR")||print "cannot open file $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;

system("rm -f $dir/count_matrix_anyBFDR_le0.05");
open(fw,">>$dir/count_matrix_anyBFDR_le0.05");
print fw "$data[0]\n";

foreach my $line(@data)
{
 if($line!~/^Prot_acc/)
 {
  @splitted_line=split("\t",$line);
  if((($splitted_line[3]<=0.05)&&($splitted_line[3] ne "NA"))||(($splitted_line[5]<=0.05)&&($splitted_line[5] ne "NA"))||(($splitted_line[7]<=0.05)&&($splitted_line[7] ne "NA")))
  {
   print fw "$line\n";
  }
 }
}
