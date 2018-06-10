#!/usr/bin/perl
use strict;

my @train;my $train;
print "Enter motif type\t";my $type=<STDIN>;chomp($type);$type=uc($type);
if($type eq "PXIXIT")
{
 $type="PxIxIT";
 @train=qw(P49790 O14511 Q9Y3S1 O95180 Q9H2S1 A1L4H1 Q6NUN7 Q9Y4F9 Q16799 Q0VAK6 Q9H252 Q14872 Q96GX8 P21333 Q9Y6V0 Q86Y37 P06731 Q8N3X6 Q9UMQ6 Q92502 P24588 Q9Y6J0 O95180 Q05193 Q9Y6M7 Q7Z418 O95644 Q13469 Q12968 Q14934 P19634 P53805);
 $train=@train;
}

if($type eq "LXVP")
{
 $type="LxVP";
 @train = qw(Q13469 Q6VAB6 Q5T0W9 Q8NFU7 O14490 P0C7V6 Q14938 Q9H4B6 P46531 Q01850 O00763 Q6IQ23 Q9NPH0 Q9HBY8 P31629 Q9HAW4 Q15842 Q6NXP0 Q7Z2Y8 O94875 Q96EV8 Q86UR5 Q9H7P9 Q9HCH0 P26368 O95180 O00429 P13861 O95644 Q12968 Q14934 Q7Z418);
 $train=@train;
}

print "Enter the entire path to the accfile to check motifs predicted in those proteins\n";my $infile=<STDIN>;chomp($infile);
open(fr,$infile)||print "cannot open $infile $!\n";#open(fr,"/data/nikhil/Desktop/acclists")||print "cannot open custom list $!\n";
my @list=<fr>;chomp(@list);

print "Enter IUPRED cut-off (motifs w IU>= will be counted)\t";my $iu=<STDIN>;chomp($iu);
print "Enter the pval cut-off (motifs w p<= will be counted)\t";my $pval=<STDIN>;chomp($pval);

my $dir="./";

print "Enter the filename w predicted motifs (one of \"top5K_hits\"/\"accessible_SLiMs_from_top5K_each\"/\"IU_le0.3\"/\"irresp_o_conservation/acc_SLiMs_pval_le0.0001\")\n";my $filename=<STDIN>;chomp($filename);
open(fr2,"$dir/$filename")||print "cannot open $filename $!\n";
my @data=<fr2>;chomp(@data);

my @splitted_line;my %all;my %prots;my $str;my %motifs;my %common;my $nomatch=0;my $flag=0;

foreach my $line(@list)
{
 $line=~s/\n//g;$line=~s/\r//g;
 $prots{$line}="";
}
my $prots=keys %prots;print "custom list has $prots proteins\n";

foreach my $line(@data)
{
 if($line!~/^pssm_p/)
 {
  @splitted_line=split("\t",$line);#print "$splitted_line[38]\n";exit();
  if(($splitted_line[12] eq $type)&&($splitted_line[38]>=$iu)&&($splitted_line[0]<=$pval))
  {
   $str=$splitted_line[5]."_".$splitted_line[8]."_".$splitted_line[9]."_".$splitted_line[10]."_".$splitted_line[12];
   $motifs{$str}="";
   $all{$splitted_line[5]}="";
  }
 }
}
my $prots=keys %all;
my $motifs=keys %motifs;print "$motifs predicted $type in $prots proteins satisfy you\n";
undef %motifs;undef %all;

foreach my $acc(sort keys %prots)
{
 foreach my $line(@data)
 {
  if($line!~/^pssm_p/)
  {
   @splitted_line=split("\t",$line);#print "$splitted_line[38]\n";exit();
   if(($acc eq $splitted_line[5])&&($splitted_line[12] eq $type)&&($splitted_line[38]>=$iu)&&($splitted_line[0]<=$pval))
   {
    $str=$splitted_line[5]."_".$splitted_line[8]."_".$splitted_line[9]."_".$splitted_line[10]."_".$splitted_line[12];
    $motifs{$str}="";
    $all{$splitted_line[5]}="";
   }
  }
 }
} 
my $motifs=keys %motifs;foreach my $motif(sort keys %motifs){print "$motif\n"}
my $all=keys %all;
print "A total of $motifs $type predicted in $all of the total $prots proteins in your input list\n";
undef %motifs;undef %all;
