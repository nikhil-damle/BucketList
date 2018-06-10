#!/usr/bin/perl
use strict;

open(fr,"PxIxIT.pssm.json")||print "cannot open pssm_PxIxIT.json $!\n";
#open(fr,"LxVP.pssm.json")||print "cannot open pssm_LxVP.json $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;my %res;my %pos;my %mat;my $res;my $pos;my $str;my @splitted_str;my $score;my $i=1;my $j=1;my @pwm;my $k;

#@splitted_line=split(/\}\,/,$data[0]);#print "$splitted_line[0]\n";
@splitted_line=split(/\]\,/,$data[0]);#print "$splitted_line[0]\n";

foreach my $ele(@splitted_line)
{
 $ele=~/\"([A-Z])\"\:/;
 $res=$1;
 $res{$res}="";#print "Dealing with res $res\n";exit();
 $str=$';
 #$str=~s/\{//g;
 $str=~s/\[//g;#print "$str\n";exit();
 @splitted_str=split(/\,/,$str);

 #foreach my $ele2(@splitted_str)
 #{
 # $ele2=~/\"(\d+)\"\:(.*)/;
 # $pos=$1;$score=$2;$score=sprintf("%.3f",$score);
 # $pos{$pos}="";
 # $mat{$res}{$pos}=$score;
 #}
 
 for($k=0;$k<scalar @splitted_str;$k++)
 {
  $pos=$k+1;$score=$splitted_str[$k];$score=sprintf("%.3f",$score);
  $pos{$pos}="";
  $mat{$res}{$pos}=$score;
 }
}
my $res=keys %res;print "There are $res residues\n";
my $pos=keys %pos;print "There are $pos positions\n";#exit();

foreach my $pos(sort {$a<=>$b} keys %pos)
{
 $pwm[0][$j]=$pos;
 $j++;
}
foreach my $res(sort keys %res)
{
 $pwm[$i][0]=$res;
 $i++;
}

foreach my $residue(sort keys %mat)
{
 foreach my $position(keys %{$mat{$residue}})
 {
  #print "$residue\t$position\t$mat{$residue}{$position}\n";
  for($i=0;$i<=$res;$i++)
  {
   if($pwm[$i][0] eq $residue)
   {
    for($j=0;$j<=$pos;$j++)
    {
     if($pwm[0][$j] eq $position)
     {
      $pwm[$i][$j]=$mat{$residue}{$position};
     }
    }
   }
  }
 }
}
for($i=0;$i<=$res;$i++)
{
 for($j=0;$j<=$pos;$j++)
 {
  print "$pwm[$i][$j]\t";
 }
 print "\n";
}
