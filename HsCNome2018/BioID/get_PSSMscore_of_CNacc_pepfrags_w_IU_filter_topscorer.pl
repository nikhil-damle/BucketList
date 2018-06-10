#!/usr/bin/perl
use strict;

print "Which motif?\t";my $motif=<STDIN>;chomp($motif);$motif=uc($motif);
print "Enter the IUPRED cut-off\t";my $cut=<STDIN>;chomp($cut);

print "Enter the dir containing acc list\n";my $dir=<STDIN>;chomp($dir);
print "Enter the name of the acc list file\t";my $inpath=<STDIN>;chomp($inpath);
print "Enter a codename to prefix the output file (\"XXX\"_CNacc_$motif\_shortIU_me$cut)\t";my $prefix=<STDIN>;chomp($prefix);

open(fr,"$dir/$inpath")||print "cannot open $dir/$inpath $!\n";
my @data=<fr>;chomp(@data);
open(fw,">>$dir/$prefix\_CNacc_$motif\_shortIU_me$cut");
my @splitted_line;my %prots;my %omit;my $i;my $acc;my $seq="";my @splitted_entry;my %scores;my %disorder;my @peporder;my $start;my $end;my @regions;my $regions;my $pass=0;my $nomatch=0;my $undef=0;my $gotit=0;my $nfac;my $cfac;

if($motif eq "PXIXIT"){$nfac=9;$cfac=2;}
if($motif eq "LXVP"){$nfac=7;$cfac=0;}

foreach my $line(@data)
{
 $line=~s/\r//;$line=~s/\n//g;
 $prots{$line}="";
}
my $prots=keys %prots;print "I got $prots to process\n";

open(fr,"Hs_proteome_usually_inaccessible_regions_Qmark_removed")||print "cannot open Hs_proteome_usually_CnInacc reg file $!\n";
my @data=<fr>;chomp(@data);

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);
 for($i=1;$i<scalar @splitted_line;$i++)
 {
  $omit{$splitted_line[0]}{$splitted_line[$i]}="";
 }
}
my $omit=keys %omit;print "Omit regions from $omit proteins\n";

local $/="\n>";
open(fr2,"Hs_rev_can_n_isoforms.fasta")||print "cannot open fasta file $!\n";
my @proteome=<fr2>;chomp(@proteome);my $len=@proteome;print "I got $len fasta seq\n";

foreach my $entry(@proteome)
{
 @splitted_line=split("\n",$entry);
 foreach my $line(@splitted_line)
 {
  if($line=~/.{2}\|(.*)\|/){$acc=$1;}
  else
  {
   $line=~s/\n//g;
   $seq=$seq.$line;
  }
 }#print "$acc\n$seq\n";
 if($acc!~/\-/)
 {
  foreach my $accin(sort keys %prots)
  {
   if($acc eq $accin)
   {
    @peporder=&peporder($acc,$seq,$motif); 
    %scores=&scores($acc,$seq,$motif);
    %disorder=&getIU($acc,$seq,$motif);

    foreach my $acc2(sort keys %omit)
    {
     if($acc eq $acc2)
     {#print "$acc\n$seq\n\n-------\n\n";
      @regions=keys %{$omit{$acc2}};$regions=@regions;#print "$acc has $regions regions\n";foreach my $reg(@regions){print "$reg\n";}exit();
      foreach my $pep(sort {$scores{$b}<=>$scores{$a}} keys %scores)
      {
       for($i=0;$i<scalar @peporder;$i++)
       {#print "Considering $i th peptide: $peporder[$i] from $acc\n";
        if($pep eq $peporder[$i])
        {
         #print "$ii\t$pep\t$scores{$pep}\t$disorder{$pep}\n";
         foreach my $reg(sort keys %{$omit{$acc}})
         {#print "Region $reg in $acc is to be omitted\n";
          if($reg=~/(\d+)\_(\d+)/)#;$start=$1;$end=$2;#print "\n\n$start ",$start-9,"\t$end ",($end-2),"\n$i\t";exit();
          {
           $start=$1;$end=$2;#print "START: $start\tEND: $end\n";
           if(($i<=($start-$nfac))||($i>=($end-$cfac)))
           {
            $pass++;#print "$i th peptide $peporder[$i] lies outside the boundaries of this region (pass = $pass)\n";
           }
          }
          else{$regions--;}
         }
         if($pass eq $regions)
         {#print "$i\t$pep\t$scores{$pep}\t$disorder{$pep}\n";
          if($disorder{$peporder[$i]}>=$cut)
          {
           print fw "$acc\t",length($seq),"\t",$i+1,"\t$peporder[$i]\t",$i+1+length($peporder[$i])-1,"\t$scores{$peporder[$i]}\t$disorder{$peporder[$i]}\n";
           $gotit=1;
          }
         }
         $pass=0;$regions=@regions;
        }
       }
       if($gotit eq 1){$gotit=0;last();}
      }
     }
     else
     {
      $nomatch++;
      if($nomatch eq $omit)
      {
       foreach my $pep(sort {$scores{$b}<=>$scores{$a}} keys %scores)
       {
        for($i=0;$i<scalar @peporder;$i++)
        {
         if(($pep eq $peporder[$i])&&($disorder{$peporder[$i]}>=$cut))
         {
          print fw "$acc\t",length($seq),"\t",$i+1,"\t$peporder[$i]\t",$i+1+length($peporder[$i])-1,"\t$scores{$peporder[$i]}\t$disorder{$peporder[$i]}\n";
          $gotit=1;
         }
        }
        if($gotit eq 1){$gotit=0;last();}
       }
      }
     }
    }
   }
   $gotit=0;
  }
 }
 $seq="";$nomatch=0;$undef=0;
}

sub scores()
{
 local $/="\n";
 my $acc=$_[0];my $seq=$_[1];my $motif=$_[2];#print "Dealing w $motif in $acc\n";
 if($motif eq "PXIXIT")
 {
  open(fr,"query_PxIxIT.pwm")||print "where's PSSM? $!\n";
 }
 if($motif eq "LXVP")
 {
  open(fr,"query_LxVP.pwm")||print "where's PSSM? $!\n";
 }
 my @pssm=<fr>;chomp(@pssm);my $reflen=@pssm;$reflen=$reflen-1;

 my $i;my $j;my @mat;my @splitted_line;my $len;my $input;my $limit;my $p;my @splitted_inp;my $k;my $ps=0;my %scores;

 for($i=1;$i<=$reflen;$i++)
 {
  @splitted_line=split("\t",$pssm[$i]);$len=@splitted_line;$len=$len-1;#print "PSSM is $len aa long\n";#exit();
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

 $limit=length($seq)-$len+1;

 for($p=0;$p<$limit;$p++)
 {
  $input=substr($seq,$p,$len);#print "checking $input from $acc\n";#exit();
  @splitted_inp=split("",$input);
  for($k=0;$k<scalar @splitted_inp;$k++)
  {
   for($i=1;$i<=$reflen;$i++)
   {
    if($splitted_inp[$k] eq $mat[$i][0])
    {
     $ps=$ps+$mat[$i][$k+1];#print "$mat[$i][0]\t$mat[$i][$k+1]\n";
    }
   }
  }
  $scores{$input}=$ps;#print "$input\t$ps\n";exit();
  $ps=0;
 }
 return %scores;
}

sub getIU()
{
 my $acc=$_[0];my $seq=$_[1];my $motif=$_[2];
 my $limit;my $len;my $i;my $cmd;my @data;my $resno;my @splitted_line;my $k;my @splitted_line2;my $line_counter=0;my $frag="";my $iu=0;my %disorder;

 if($motif eq "PXIXIT"){$len=12;}
 if($motif eq "LXVP"){$len=7;}

 #system("rm -f ./tmp.fasta");
 open(tmp,">tmp.fasta"); print tmp ">$acc\n$seq\n"; close tmp;

 $limit=length($seq)-$len+1;

 $cmd=`iupred ./tmp.fasta short`;
 @data=split(/\n/,$cmd);

 for($i=0;$i<$limit+9;$i++)
 {
  $data[$i]=~s/^\s+//;
  if($data[$i]!~/^#/)
  {
   @splitted_line=split(/\s+/,$data[$i]);
   $resno=$splitted_line[0];

   for($k=$i;$k<$i+$len;$k++)
   {
    $data[$k]=~s/^\s+//;
    @splitted_line2=split(/\s+/,$data[$k]);
    $line_counter++;#print "appending $splitted_line2[1]\n";
    $frag=$frag.$splitted_line2[1];
    $iu=$iu+$splitted_line2[2];
   }
   $iu=$iu/$line_counter;$iu=sprintf("%.3f",$iu);
   $disorder{$frag}=$iu;
   $line_counter=0;$iu=0;$frag="";#last();#exit();
  }
 }
 return %disorder;
}

sub peporder()
{
 my $acc=$_[0];my $seq=$_[1];my $motif=$_[2];
 my $limit;my $len;my $i;my @peps;

 if($motif eq "PXIXIT"){$len=12;}
 if($motif eq "LXVP"){$len=7;}

 $limit=length($seq)-$len+1;

 for($i=0;$i<$limit;$i++)
 {
  push(@peps,substr($seq,$i,$len));
 }
 #print "$acc\n$seq\n\n---------\n\n@peps\n\n------------------\n\n";
 return @peps;
}
