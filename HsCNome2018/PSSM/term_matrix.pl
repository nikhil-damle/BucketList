#!/usr/bin/perl
use strict;

my @splitted_line;my %pred;my $nomatch=0;my %nomatches;my %prs;my %nomatches2;my $BP;my @data;my %prots;my $acc;my %BPs;my %match;my %match2;my $gene;my $description;my $str;
my @mat;my $i;my $j;my $k;my @BPs_arr;my $genes;my @genes;my $match=0;my %commons;my @procs;my $procs;my %hash;my @splitted_str;

my $dir="/HsCNome_analyses_Feb2018/GO_analyses/BPs";
#my $dir="/HsCNome_analyses_Feb2018/GO_analyses/MFs";
#my $dir="/HsCNome_analyses_Feb2018/GO_analyses/CCs";
opendir(DIR,"$dir")||die $!;

system("rm -f $dir/gene_term_associations");
open(fw,">>$dir/gene_term_associations");

foreach my $file(readdir DIR)
{
 if($file=~/(.*)\.txt/)
 {
  $BP=$1;print "$BP\n";
  open(fr,"$dir/$file")||print "cannot open $file $!\n";
  @data=<fr>;chomp(@data);

  foreach my $line(@data)
  {
   @splitted_line=split("\t",$line);

   $acc=$splitted_line[1];
   $splitted_line[2]=~/\;(.*)\;/g;
   $gene=$1;
   $description=$splitted_line[3];
   $str=$acc."_".$gene."_".$description;

   $BPs{$BP}{$str}="";
   $prots{$acc}="";
   push(@{$hash{$str}},$BP);

   #if($line=~/UniProtKB=(.*)/)
   #{
   # $acc=$1;#print "$acc\n";exit();
   # $BPs{$BP}{$acc}="";
   # $prots{$acc}="";
   # push(@{$hash{$acc}},$BP);
   #}
  }
 }
}
my $prots=keys %prots;print "\nAforementioned BPs a/c for $prots prots\n";

foreach my $str(sort keys %hash)
{
 @splitted_str=split(/\_/,$str);$acc=$splitted_str[0];$gene=$splitted_str[1];$description=$splitted_str[2];
 $procs=@{$hash{$str}};
 print fw "$acc\t$gene\t$procs\t@{$hash{$str}}\n";
 #foreach my $line(@mapping)
 #{
  #@splitted_line=split("\t",$line);
  #if($acc eq $splitted_line[0])
  #{
   #print "$acc\t$splitted_line[1]\t$procs\t@{$hash{$acc}}\n";
  #}
 #}
}

@BPs_arr=keys %BPs;#print "@BPs\n";

for($i=0;$i<=scalar @BPs_arr;$i++)
{
 $mat[0][$i+1]=$BPs_arr[$i];
 $mat[$i+1][0]=$BPs_arr[$i];
}

for($i=0;$i<scalar @BPs_arr;$i++)
{
 for($j=0;$j<scalar @BPs_arr;$j++)
 {
  foreach my $acc1(sort keys %{$BPs{$BPs_arr[$i]}})
  {
   foreach my $acc2(sort keys %{$BPs{$BPs_arr[$j]}})
   {
    if($acc1 eq $acc2)
    {
     $match++;
     if($i ne $j)
     {
      $commons{$acc1}++;
      #foreach my $line(@mapping)
      #{
       #@splitted_line=split("\t",$line);
       #if($acc1 eq $splitted_line[0])
       #{
        ##print "$splitted_line[1]\tis common b/w $BPs_arr[$i] and $BPs_arr[$j]\n";
       #}
      #}
     }
    }
   }
  }
  $mat[$i+1][$j+1]=$match;
  $match=0;
 }
}

for($i=0;$i<scalar @mat;$i++)
{
 for($j=0;$j<scalar @mat;$j++)
 {
#  print "$mat[$i][$j]\t";
 }
# print "\n";
}
my $commons=keys %commons;print "$commons prots annotated w >1 BPs\n";

=headforeach my $acc(sort keys %commons)
{#print "$acc\n";
 foreach my $line(@mapping)
 {
  @splitted_line=split("\t",$line);
  if($acc eq $splitted_line[0])
  {
   #print "CN\tpp\t$splitted_line[1]\n";
  }
 }
}
=cut
