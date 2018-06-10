########### BASIC STATS OF THE PREDICTED PSSM HITS ############
#!/usr/bin/perl
use strict;

print "Enter the file (entire path incl the name) w predicted Cn-binding SLiMs\n";my $infile=<STDIN>;chomp($infile);

open(fr,$infile)||print "cannot open input $!\n";
my @ND_filtered=<fr>;chomp(@ND_filtered);

my @splitted_line;my %ND_filtered;my %prots;my @slims;my %pxixit_only;my %lxvp_only;my %both;my %nothing;my $str="";my %px;my %lx;

foreach my $line(@ND_filtered)
{
 if($line!~/^pssm_p/)
 {
  @splitted_line=split("\t",$line);
  $ND_filtered{$splitted_line[5]}="";
  $prots{$splitted_line[5]}{$splitted_line[12]}="";
  if($splitted_line[12] eq "PxIxIT")
  {
   $str=$splitted_line[5]."_".$splitted_line[8]."_".$splitted_line[9]."_".$splitted_line[10]."_".$splitted_line[0]."_".$splitted_line[12];
   $px{$str}="";
  }
  if($splitted_line[12] eq "LxVP")
  {
   $str=$splitted_line[5]."_".$splitted_line[8]."_".$splitted_line[9]."_".$splitted_line[10]."_".$splitted_line[0]."_".$splitted_line[12];
   $lx{$str}="";
  }
  $str="";
 }
}
my $ND_filtered=keys %ND_filtered;print "$ND_filtered prots were retained by Norman after applynig some filters\n";#exit();
my $px=keys %px;my $lx=keys %lx;print "You have $px PxIxIT SLiMs and $lx LxVP SLiMs predicted\n";

foreach my $acc(keys %prots)
{
 @slims=keys %{$prots{$acc}};#print "$slims[0]\n";exit();
 if((scalar @slims eq 1)&&($slims[0] eq "LxVP"))
 {
  $lxvp_only{$acc}="";
 }
 if((scalar @slims eq 1)&&($slims[0] eq "PxIxIT"))
 {
  $pxixit_only{$acc}="";
 }
 if((scalar @slims eq 2)&&(($slims[0] eq "PxIxIT")||($slims[0] eq "LxVP")))
 {
  $both{$acc}="";
 }
 if((scalar @slims eq 1)&&($slims[0] eq ""))
 {
  $nothing{$acc}="";
 }
}
my $Px_only=keys %pxixit_only;print "$Px_only have only PxIxIT\n";
my $Lx_only=keys %lxvp_only;print "$Lx_only have only LxVP\n";
my $both=keys %both;print "$both have both PxIxIT and LxVP\n";

foreach my $acc(sort keys %both){print "$acc\n";}
