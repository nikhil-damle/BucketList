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
my @splitted_line;my %prots;my %omit;my $i;my $acc;my $seq="";my @splitted_entry;my %scores;my %disorder;my @peporder;my $start;my $end;my @regions;my $regions;my $pass=0;my $nomatch=0;my $undef=0;

foreach my $line(@data)
{
 $line=~s/\r//;$line=~s/\n//g;
 $prots{$line}="";
}
my $prots=keys %prots;print "I got $prots to process\n";              # Total # of input proteins

open(fr,"Hs_proteome_usually_inaccessible_regions")||print "cannot open Acclist_to_b_processed $!\n";
my @data=<fr>;chomp(@data);

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);
 for($i=1;$i<scalar @splitted_line;$i++)
 {
  $omit{$splitted_line[0]}{$splitted_line[$i]}="";
 }
}
my $omit=keys %omit;print "Omit regions from $omit proteins\n";       # Proteins containing CN-inacc regions

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
   $seq=$seq.$line;                                     # acc and seq of each protein in the proteome
  }
 }
 if($acc!~/\-/)                                         # if the above acc is not that of an isoform
 {
  foreach my $accin(sort keys %prots)
  {
   if($acc eq $accin)                                   # and if it's one of the acc in the input
   {
    @peporder=&peporder($acc,$seq,$motif);              # @peporder has all 12mer or 7mers in that protein from N to C
    %scores=&scores($acc,$seq,$motif);                  # %scores{12mer or 7mer} = PSSM score (PSSM score of a given pep is independent of the N to C order)
    %disorder=&getIU($acc,$seq,$motif);

    foreach my $acc2(sort keys %omit)
    {
     if($acc eq $acc2)                                  # if that acc contains a Cn-inacc region
     {
      @regions=keys %{$omit{$acc2}};$regions=@regions;#print "$acc has $regions regions\n";foreach my $reg(@regions){print "$reg\n";}exit();
      for($i=0;$i<scalar @peporder;$i++)                # Reading the  protein from N to C
      {#print "Considering $i th peptide: $peporder[$i] from $acc\n";
       foreach my $reg(sort keys %{$omit{$acc}})
       {
        #print "Region $reg in $acc is to be omitted\n";
        if($reg=~/(\d+)\_(\d+)/)
        {
         $start=$1;$end=$2;#print "START: $start\tEND: $end\n";
         if(($i<=($start-9))||($i>=($end-2)))           # Ensuring that the frag in question DOES NOT overlap considerably w the CN-inacc reg
         {
          $pass++;                                      # counter pass increases only if the pep lies outside the Cn-inacc region
          #print "$i th peptide $peporder[$i] lies outside the boundaries of this region (pass = $pass)\n";
         }
        }
        elsif($reg=~/\?/g){$regions--;}#print "Total regions to be scanned remain: $regions\n";}   # $region decreases by 1 with each CN-inacc region checked
       }
       if($pass eq $regions)
       {
        #print "$disorder{$peporder[$i]}\t";exit();
        if($disorder{$peporder[$i]}>=$cut)               # now if this pep lying outside of the CN-inacc reg passes the IU cut-off
        {
         #print "$acc\t$peporder[$i]\t$scores{$peporder[$i]}\t$disorder{$peporder[$i]}\t";exit();
         print fw "$acc\t",length($seq),"\t",$i+1,"\t$peporder[$i]\t",$i+1+length($peporder[$i])-1,"\t$scores{$peporder[$i]}\t$disorder{$peporder[$i]}\n";
        }
       }
       $pass=0;$regions=@regions;
      }
      #print "\n";#exit();
     }
     else
     {
      $nomatch++;
      if($nomatch eq $omit)                              # if the protein contains no CN-inacc regions simply get the top-scoring frag that fits IU cut-off
      {
       for($i=0;$i<scalar @peporder;$i++)
       {
        if($disorder{$peporder[$i]}>=$cut)
        {
         print fw "$acc\t",length($seq),"\t",$i+1,"\t$peporder[$i]\t",$i+1+length($peporder[$i])-1,"\t$scores{$peporder[$i]}\t$disorder{$peporder[$i]}\n";
        }
       }
      }
     }
    }
   }
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

 for($i=0;$i<$limit+9;$i++)                    # first 9 lines of the output start w # and need not be read
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
