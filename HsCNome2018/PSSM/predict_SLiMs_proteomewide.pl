#!/usr/bin/perl
use strict;

my @splitted_entry;my $seq="";my $limit;my $len;my $acc;my $inpath;my @data;my $p;my %scores;
my $i;my $j;my @mat;my @splitted_line;my $input;my $k;my @splitted_inp;my $ps=0;

open(fw,">>./PSSM_scores");

print "Enter the motif type(PxIxIT/LxVP):\t";my $type=<STDIN>;chomp($type);
if(($type eq "PxIxIT")||($type eq "pxixit")){$type="PxIxIT";}
if(($type eq "LxVP")||($type eq "lxvp")){$type="LxVP";}
print "Enter the score cut-off:\t";my $cut=<STDIN>;chomp($cut);

open(fr,"./query_$type\.pwm")||print "where's PSSM? $!\n";
my @pssm=<fr>;chomp(@pssm);my $reflen=@pssm;$reflen=$reflen-1;

######## Reading PWMs ########

for($i=1;$i<=$reflen;$i++)
{
 @splitted_line=split("\t",$pssm[$i]);$len=@splitted_line;$len=$len-1;#print "$len\n";
 $mat[$i][0]=$splitted_line[0];
}
for($j=1;$j<=$len;$j++)
{
 $mat[0][$j]="pos".$j;
}
for($i=1;$i<=$reflen;$i++)
{
 @splitted_line=split("\t",$pssm[$i]);
 for($j=1;$j<=$len;$j++)
 {
  $mat[$i][$j]=$splitted_line[$j];
 }
}
#for($i=0;$i<=$reflen;$i++)
#{
# for($j=0;$j<=$len;$j++)
# {
#  print "$mat[$i][$j]\t";
# }
# print "\n";
#}
#exit();
########## Reading PWMs done #######

print "Do u have a single stretch (press 1) or a file (press 2) ?\t";my $choice=<STDIN>;chomp($choice);
if(($choice!~/[0-9]/)||(($choice ne 1)&&($choice ne 2))){print "Don't be mischievous; press 1 or 2\n";}

if($choice eq 1)
{
 print "Enter your query aa stretch\t";$input=<STDIN>;chomp($input);
 $input=uc($input);$input=~s/\n//;$input=~s/\r//;
 if($input=~/[BJOUXZ]/){print "Input seq has atypical amino acids\n";exit();}
 if((length($input)<$len)||(length($input)>$len)){print "Give me a $len mer stretch\n";exit();}
 @splitted_inp=split("",$input);
 for($k=0;$k<scalar @splitted_inp;$k++)
 {
  for($i=1;$i<=$reflen;$i++)
  {
   if($splitted_inp[$k] eq $mat[$i][0])
   {
    #print "$splitted_inp[$k] in row $i @ pos $k+1 i.e. $mat[$i][$k+1]\n";
    $ps=$ps+$mat[$i][$k+1];
   }
  }
 }
 print "Combined score for $input is $ps\n";
 print fw "$input\t$ps\n";
}
if($choice eq 2)
{
 print "Enter the entire path to fasta file\n";$inpath=<STDIN>;chomp($inpath);
 local $/="\n>";
 open(fr2,$inpath)||print "Where's $inpath $!\n";
 @data=<fr2>;chomp(@data);print "You gave me ", scalar @data, " seq\n";#exit();

 foreach my $entry(@data)
 {
  $entry=~s/>//;
  @splitted_entry=split("\n",$entry);
  $splitted_entry[0]=~/.{2}\|(.*)\|/;$acc=$1;print "Dealing with $acc\n";   # Modify here according to your file's header
  for($i=1;$i<scalar @splitted_entry;$i++)
  {
   $seq=$seq.$splitted_entry[$i];
  }
  $limit=length($seq)-$len+1;print "You gave me ", length($seq)," aa long seq; you will have $limit $len long $type stretches in it\n";#exit();
  for($p=0;$p<$limit;$p++)
  {
   $input=substr($seq,$p,$len);#print "I have $input which is ",length($input)," aa long & starts @ ",$p+1," in $acc\n";
   $input=uc($input);$input=~s/\n//;$input=~s/\r//;

   if($input=~/[BJOUXZ]/){print "Input seq has atypical amino acids\n";exit();}
   if((length($input)<$len)||(length($input)>$len)){print "Give me a $len mer stretch\n";exit();}

   @splitted_inp=split("",$input);
   for($k=0;$k<scalar @splitted_inp;$k++)
   {
    for($i=1;$i<=$reflen;$i++)
    {
     if($splitted_inp[$k] eq $mat[$i][0])
     {
      #print "$splitted_inp[$k] in row $i @ pos $k+1 i.e. $mat[$i][$k+1]\n";
      $ps=$ps+$mat[$i][$k+1];
     }
    }
   }
   #print "Combined score for $input is $ps\n";
   $scores{$input}=$ps;
   $ps=0;
  }
  foreach my $pep(sort {$scores{$b}<=>$scores{$a}} keys %scores)
  {
   if($scores{$pep} >= $cut)
   {
    print fw "$acc\t$pep\t$scores{$pep}\n";#last();
   }
  }
  $seq="";undef %scores;#exit();
 }
}
