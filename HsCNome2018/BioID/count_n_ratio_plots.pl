#!/usr/bin/perl
use strict;

print "Enter the working dir\t";my $dir=<STDIN>;chomp($dir);

#my $dir="BioID";
#my $dir="SWATH";
#my $dir="DDA";

open(fr,"$dir/count_matrix_anyBFDR_le0.05_w_mapping")||print "cannot open cout file $!\n";
my @data=<fr>;chomp(@data);

system("rm -f $dir/count_matrix_anyBFDR_le0.05_w_mapping_for_plotting");
open(fw,">>$dir/count_matrix_anyBFDR_le0.05_w_mapping_for_plotting");

print fw "GeneID\tMapped_acc\tGene\tWT_count\tWT_BFDR\tNIR_count\tNIR_BFDR\tWF_count\tWF_BFDR\tWT_NIR_ratio\tWT_WF_ratio\n";

my @splitted_line;my %hash;my $i;my $wt_nir;my $wt_wf;

for($i=3;$i<scalar @data;$i++)
{
 @splitted_line=split("\t",$data[$i]);
 if($splitted_line[1] ne "Not mapped")
 {
  if($splitted_line[5] eq "0")
  {
   $wt_nir=$splitted_line[3];
  }
  if($splitted_line[5] ne "0")
  {
   $wt_nir=$splitted_line[3]/$splitted_line[5];
   $wt_nir=sprintf("%.3f",$wt_nir);
  }
  if($splitted_line[7] eq "0")
  {
   $wt_wf=$splitted_line[3];
  }
  if($splitted_line[7] ne "0")
  {
   $wt_wf=$splitted_line[3]/$splitted_line[7];
   $wt_wf=sprintf("%.3f",$wt_wf);
  }
  print fw "$data[$i]$wt_nir\t$wt_wf\n";#exit();
 }
}
close fw;

open(fr,"$dir/count_matrix_anyBFDR_le0.05_w_mapping_for_plotting")||print "cannot open the file $!\n";
my @data=<fr>;chomp(@data);

getlog(\@data,"$dir/count_matrix_anyBFDR_le0.05_w_mapping_for_plotting_w_log2");

sub getlog()
{
 my @data=@{$_[0]};my @splitted_line;my $file=$_[1];my $log2nir;my $log2wf;my $i;
 system("rm -f $file");
 open(fw,">>$file")||print "cannot open $file $!\n";
 print fw "$data[0]\tWT_NIR_log2\tWT_WF_log2\n";
 for($i=1;$i<scalar @data;$i++)
 {
  @splitted_line=split("\t",$data[$i]);
  if(($splitted_line[9] ne "0")&&($splitted_line[9] ne "0.000")){$log2nir=log($splitted_line[9])/log(2);}
  else{$log2nir=0;}
  if(($splitted_line[10] ne "0")&&($splitted_line[10] ne "0.000")){$log2wf=log($splitted_line[10])/log(2);}
  else{$log2wf=0;}#print "$splitted_line[9]\t$splitted_line[10]\t$log2wf\n";}
  print fw "$data[$i]\t$log2nir\t$log2wf\n";
 }
}
