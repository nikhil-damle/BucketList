#!/usr/bin/perl
use strict;

my $dir1="SWATH";
my $dir2="DDA";

my @splitted_line;my %swath;my %dda;my %signal;my %both;my %swath_signal;my %dda_signal;my @acc;

print "Enter the log2(WT:Mut) ratio cut-off (enter na for not-cut-off)\t";my $cut=<STDIN>;chomp($cut);if($cut eq "na"){$cut=-999;}

open(fr,"$dir1/count_matrix_anyBFDR_le0.05_w_mapping_for_plotting_w_log2")||print "cannot open count_matrix_anyBFDR_le0.05_w_mapping_for_plotting_w_log2 $!\n";
my @data=<fr>;chomp(@data);

system("rm -f SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2");
system("rm -f SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2_ratio_gt$cut");

open(fw,">>SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2");
open(fw2,">>SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2_ratio_gt$cut");

open(fw,">>$dir1/count_matrix_WT_BFDR_le0.05_w_mapping_for_plotting_w_log2");

foreach my $line(@data)
{
 if($line!~/^GeneID/)
 {
  @splitted_line=split("\t",$line);
  if(($splitted_line[4]<=0.05)&&($splitted_line[4] ne "NA"))
  {
   print fw "$line\tSWATH\n";
   if($splitted_line[1]=~/\,/)
   {
    $splitted_line[1]=~s/\s+//g;
    @acc=split(/\,/,$splitted_line[1]);#print "@acc\n";exit();
    foreach my $acc(sort @acc)
    {
     $swath{$acc}="";
     $both{$acc}="";
     if(($splitted_line[11] > $cut)||($splitted_line[12]>$cut))
     {
      print fw2 "$line\tSWATH\n";
      $swath_signal{$acc}="";
      $signal{$acc}="";
     }
    }
   }
   else
   {
    $swath{$splitted_line[1]}="";
    $both{$splitted_line[1]}="";
    if(($splitted_line[11] > $cut)||($splitted_line[12]>$cut))
    {
     print fw2 "$line\tSWATH\n";
     $swath_signal{$splitted_line[1]}="";
     $signal{$splitted_line[1]}="";
    }
   }
  }
 }
}
my $swath=keys %swath;
my $swath_signal=keys %swath_signal;
print "$swath proteins in SWATH expt had WT-BFDR<=0.05 of which $swath_signal proteins showed reduced BioID signal with either of the CN-muts\n";

open(fr,"$dir2/count_matrix_anyBFDR_le0.05_w_mapping_for_plotting_w_log2")||print "cannot open count_matrix_anyBFDR_le0.05_w_mapping_for_plotting_w_log2 $!\n";
my @data=<fr>;chomp(@data);

open(fw,">>$dir2/count_matrix_WT_BFDR_le0.05_w_mapping_for_plotting_w_log2");

foreach my $line(@data)
{
 if($line!~/^GeneID/)
 {
  @splitted_line=split("\t",$line);
  if(($splitted_line[4]<=0.05)&&($splitted_line[4] ne "NA"))
  {
   print fw "$line\tDDA\n";
   if($splitted_line[1]=~/\,/)
   {
    $splitted_line[1]=~s/\s+//g;
    @acc=split(/\,/,$splitted_line[1]);#print "@acc\n";exit();
    foreach my $acc(sort @acc)
    {
     $dda{$acc}="";
     $both{$acc}="";
     if(($splitted_line[11] > $cut)||($splitted_line[12]>$cut))
     {
      print fw2 "$line\tDDA\n";
      $dda_signal{$acc}="";
      $signal{$acc}="";
     }
    }
   }
   else
   {
    $dda{$splitted_line[1]}="";
    $both{$splitted_line[1]}="";
    if(($splitted_line[11] > $cut)||($splitted_line[12]>$cut))
    {
     print fw2 "$line\tDDA\n";
     $dda_signal{$splitted_line[1]}="";
     $signal{$splitted_line[1]}="";
    }
   }
  }
 }
}
my $dda=keys %dda;
my $dda_signal=keys %dda_signal;
print "$dda proteins in DDA expt had WT-BFDR<=0.05 of which $dda_signal proteins showed reduced BioID signal with either of the CN-muts\n";

my $both=keys %both;
my $signal=keys %signal;
print "$both proteins in SWATH+DDA expts had WT-BFDR<=0.05 of which $signal proteins showed reduced BioID signal with either of the CN-muts in either of the SWATH or DDA expt\n";

#foreach my $acc(sort keys %both)
#foreach my $acc(sort keys %signal)
#{
# print "$acc\n";
#}
