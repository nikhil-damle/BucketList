#!/usr/bin/perl
use strict;

local $/="//\n";
open(fr,"Hs_proteome.txt")||print "cannot open Hs_proteome.txt $!\n";
my @data=<fr>;chomp(@data);my $prots=@data;print "Uniprot on Nov2017 has $prots reviewed proteins\n";

system("rm -f Hs_proteome_usually_inaccessible_regions Acclist_to_b_processed");
open(fw,">>Hs_proteome_usually_inaccessible_regions");
open(fw2,">>Acclist_to_b_processed");

my @splitted_data;my @splitted_line;my @accarr;my $acc;my %omit;my $omit;my @peps;my $flag=0;my %secreted;my %process;

foreach my $entry(@data)
{
 @splitted_data=split("\n",$entry);
 foreach my $line(@splitted_data)
 {
  if($line=~/^AC   /)
  {
   @splitted_line=split(/\;/,$');
   foreach my $acc(@splitted_line)
   {
    push(@accarr,$acc);
   }
   $acc=$accarr[0];
  }
  if($line=~/^KW.*Secreted/)
  {
   $flag=1;
  } 
 }
 $acc=$accarr[0];
 if($flag eq 1){$secreted{$acc}="";}
 if($flag eq 0)
 {
  foreach my $line(@splitted_data)
  {
   if($line=~/^FT   /)
   {
    @splitted_line=split(/\s+/,$line);
    if(($splitted_line[1] ne "")&&($splitted_line[2] ne "")&&($splitted_line[3] ne "")&&(($splitted_line[4]=~/Extracellular/i)||($splitted_line[4]=~/Mitochondrion/i)||($splitted_line[4]=~/Vesicular/i)||($splitted_line[4]=~/Lumenal/i)||($splitted_line[4]=~/Mitochondrial intermembrane/i)||($splitted_line[4]=~/Lysosome/i)||($splitted_line[4]=~/Endosome/i)||($splitted_line[1] eq "TRANSMEM")||($splitted_line[1]=~/INTRAMEM/i)||($splitted_line[4]=~/Peroxisome/i)||($splitted_line[4]=~/Perinuclear space/i)||($splitted_line[4]=~/Intragranular/i)||($splitted_line[4]=~/Vacuolar/i)||($splitted_line[4]=~/Mitochondrial matrix/i)||($splitted_line[4]=~/Exoplasmic loop/i)))
    {
     $omit=$splitted_line[2]."_".$splitted_line[3];
     $omit{$acc}{$omit}="";
    }
    else
    {
     $process{$acc}="";
    }
   }
  }
 }
 undef @accarr;$flag=0;#last();
}
my $secreted=keys %secreted;print "$secreted proteins are annotated as secreted\n";
foreach my $acc(sort keys %omit)
{
 print fw "$acc\t";
 foreach my $limits(sort keys %{$omit{$acc}})
 {
  print fw "$limits\t";
 }
 print fw "\n";
}
foreach my $acc(sort keys %process)
{
 print fw2 "$acc\n";
}
