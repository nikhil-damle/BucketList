#!/usr/bin/perl
use strict;

print "Enter the path to working directory\n";my $dir=<STDIN>;chomp($dir);
print "Enter ProP-PD fasta file (just file name, you already gave me its location above)\t";my $file=<STDIN>;chomp($file);
print "Your o/p redirected to $dir/motif_matrix\n";
open(fr,"$dir/$file")||print "cannot open all_counts_12Oct_corrected_with_pos.fasta $!\n";
my @data=<fr>;chomp(@data);
system("rm -f $dir/motif_matrix");
open(fw,">>$dir/motif_matrix");
my @splitted_line;my %peps;my $acc;my $k=1;
my $m1flag=0;my $m2flag=0;my $m3flag=0;my $m4flag=0;my $m5flag=0;my $m6flag=0;my $m7flag=0;my $m8flag=0;my $m9flag=0;my $m10flag=0;my $m11flag=0;my $m12flag=0;my $m13flag=0;my $m14flag=0;my $m15flag=0;my $m16flag=0;my $m17flag=0;my $m18flag=0;my $m19flag=0;my $m20flag=0;my $m21flag=0;my $m22flag=0;my $m23flag=0;

foreach my $line(@data)
{
 if($line=~/^>(.*)\|.*/)
 {
  $acc=$1;
 }
 else
 {
  #if($line=~/[PI][A-Z][ILVF][A-Z][A-Z][A-Z]/)
  if($line=~/[PI][^PG][ILVF][^PG][ILVF][TSHDEQNKR]/)             # No variation
  {
   $line=~s/\n//g;
   $m1flag=1;
  }
  #if($line=~/[PI][A-Z][ILVF][A-Z][A-Z][TSHDEQNKRILVF]/)
  if($line=~/[^PI][^PG][ILVF][^PG][ILVF][TSHDEQNKR]/)            # Pos1
  {
   $line=~s/\n//g;
   $m2flag=1;
  }
  #if($line=~/[PI][A-Z][ILVF][A-Z][A-Z][TSHDEQNKR]/)
  if($line=~/[PI][PG][ILVF][^PG][ILVF][TSHDEQNKR]/)              # Pos2
  {
   $line=~s/\n//g;
   $m3flag=1;
  }
  #if($line=~/[PI][^PG][ILVF][^PG][A-Z][TSHDEQNKR]/)
  if($line=~/[PI][^PG][^ILVF][^PG][ILVF][TSHDEQNKR]/)            # Pos3
  {
   $line=~s/\n//g;
   $m4flag=1;
  }
  #if($line=~/[A-Z][A-Z][ILVF][A-Z][ILVF][TSHDEQNKR]/)
  if($line=~/[PI][^PG][ILVF][PG][ILVF][TSHDEQNKR]/)              # Pos4
  {
   $line=~s/\n//g;
   $m5flag=1;
  }
  #if($line=~/[PI][A-Z][A-Z][A-Z][ILVF][TSHDEQNKR]/)
  if($line=~/[PI][^PG][ILVF][^PG][^ILVF][TSHDEQNKR]/)            # Pos5
  {
   $line=~s/\n//g;
   $m6flag=1;
  }
  #if($line=~/[PI][^PG][A-Z][^PG][ILVF][TSHDEQNKR]/)
  if($line=~/[PI][^PG][ILVF][^PG][ILVF][^TSHDEQNKR]/)            # Pos6
  {
   $line=~s/\n//g;
   $m7flag=1;
  }
  #if($line=~/[PI][^PG][A-Z][^PG][A-Z][ILVF]/)
#  if($line=~/[PI][^PG][^ILVF][^PG][^ILVF][TSHDEQNKR]/)           # Pos3,5
#  {
#   $line=~s/\n//g;
#   $m8flag=1;
#  }
  #if($line=~/[A-Z][^PG][ILVF][^PG][ILVF][TSHDEQNKR]/)
#  if($line=~/[PI][^PG][ILVF][^PG][^ILVF][^TSHDEQNKR]/)           # Pos5,6
#  {
#   $line=~s/\n//g;
#   $m9flag=1;
#  }
  #if($line=~/[PI][^PG][ILVF][^PG][ILVF][A-Z]/)
#  if($line=~/[PI][^PG][^ILVF][^PG][ILVF][^TSHDEQNKR]/)           # Pos3,6
#  {
#   $line=~s/\n//g;
#   $m10flag=1;
#  }
  #if($line=~/[PI][A-Z][ILVF][A-Z][ILVF][TSHDEQNKR]/)
#  if($line=~/[^PI][^PG][^ILVF][^PG][ILVF][TSHDEQNKR]/)           # Pos1,3
#  {
#   $line=~s/\n//g;
#   $m11flag=1;
#  }
  #if($line=~/[PI][A-Z][ILVF][^PG][ILVF][TSHDEQNKR]/)
  if($line=~/L[A-Z]VP/)                                          # LxVP
  {
   $line=~s/\n//g;
   $m12flag=1;
  }
  if($line=~/P[A-Z]{3}[IV][TDH]/)                                # Peti's Px1
  {
   $line=~s/\n//g;
   $m13flag=1;
  }
  if($line=~/P[A-Z]{4}[IV][TDH]/)                                # Peti's Px2
  {
   $line=~s/\n//g;
   $m14flag=1;
  }
  if($line=~/[NQDESTRH][YTDILVF]L[A-Z]VP/)                        # Peti's strict-Lx
  {
   $line=~s/\n//g;
   $m15flag=1;
  }
  if($line=~/[NQDESTRH][YTDILVF]L[A-Z][VPL][PK]/)                 # Peti's extended-Lx
  {
   $line=~s/\n//g;
   $m16flag=1;
  }
=head  if($line=~/[PI][A-Z][ILVF][A-Z][A-Z][ILVF]/)
  {
   $line=~s/\n//g;
   $m13flag=1;
  }
  if($line=~/[PI][^PG][ILVF][A-Z][ILVF][TSHDEQNKR]/)
  {
   $line=~s/\n//g;
   $m14flag=1;
  }
  if($line=~/[PI][^PG][ILVF][^PG][ILVF][TSHDEQNKR]/)
  {
   $line=~s/\n//g;
   $m15flag=1;
  }
  if($line=~/P[^PG][ILVF][^PG][ILVF][TSHDEQNKR]/)
  {
   $line=~s/\n//g;
   $m16flag=1;
  }
  if($line=~/[A-Z][A-Z][ILVF][A-Z][ILVF][ILVF]/)
  {
   $line=~s/\n//g;
   $m17flag=1;
  }
  if($line=~/[PI][A-Z][A-Z][A-Z][ILVF][ILVF]/)
  {
   $line=~s/\n//g;
   $m18flag=1;
  }
  if($line=~/[PI][^PG][A-Z][^PG][A-Z][TSHDEQNKR]/)
  {
   $line=~s/\n//g;
   $m19flag=1;
  }
  if($line=~/[PI][A-Z][ILVF][A-Z][ILVF][ILVF]/)
  {
   $line=~s/\n//g;
   $m20flag=1;
  }
  if($line=~/[PI][^PG][ILVF][^PG][ILVF][ILVF]/)
  {
   $line=~s/\n//g;
   $m21flag=1;
  }
  if($line=~/I[^PG][ILVF][^PG][ILVF][TSHDEQNKR]/)
  {
   $line=~s/\n//g;
   $m22flag=1;
  }
  if($line=~/[A-Z][ILVF][A-Z][ILVFP]P/)
  {
   $line=~s/\n//g;
   $m23flag=1;
=cut}
  #if(($m1flag ne 0)||($m2flag ne 0)||($m3flag ne 0)||($m4flag ne 0)||($m5flag ne 0)||($m6flag ne 0)||($m7flag ne 0)||($m8flag ne 0)||($m9flag ne 0)||($m10flag ne 0)||($m11flag ne 0)||($m12flag ne 0)||($m13flag ne 0)||($m14flag ne 0)||($m15flag ne 0)||($m16flag ne 0)||($m17flag ne 0)||($m18flag ne 0)||($m19flag ne 0)||($m20flag ne 0)||($m21flag ne 0)||($m22flag ne 0)||($m23flag ne 0))
  #if(($m1flag ne 0)||($m2flag ne 0)||($m3flag ne 0)||($m4flag ne 0)||($m5flag ne 0)||($m6flag ne 0)||($m7flag ne 0)||($m8flag ne 0)||($m9flag ne 0)||($m10flag ne 0)||($m11flag ne 0)||($m12flag ne 0)||($m13flag ne 0)||($m14flag ne 0))#||($m15flag ne 0)||($m16flag ne 0)||($m17flag ne 0)||($m18flag ne 0)||($m19flag ne 0)||($m20flag ne 0)||($m21flag ne 0)||($m22flag ne 0)||($m23flag ne 0))
  #{
   #print fw "$acc\t$line\tP$k\t$m1flag\t$m2flag\t$m3flag\t$m4flag\t$m5flag\t$m6flag\t$m7flag\n";#\t$m8flag\t$m9flag\t$m10flag\t$m11flag\t$m12flag\t$m13flag\t$m14flag\t$m15flag\t$m16flag\n";#$m17flag\t$m18flag\t$m19flag\t$m20flag\t$m21flag\t$m22flag\t$m23flag\n";
   print fw "$acc\t$line\tP$k\t$m1flag\t$m2flag\t$m3flag\t$m4flag\t$m5flag\t$m6flag\t$m7flag\t$m12flag\t$m13flag\t$m14flag\t$m15flag\t$m16flag\n";#$m17flag\t$m18flag\t$m19flag\t$m20flag\t$m21flag\t$m22flag\t$m23flag\n";
  #}
  $m1flag=0;$m2flag=0;$m3flag=0;$m4flag=0;$m5flag=0;$m6flag=0;$m7flag=0;$m8flag=0;$m9flag=0;$m10flag=0;$m11flag=0;$m12flag=0;$m13flag=0;$m14flag=0;$m15flag=0;$m16flag=0;#$m17flag=0;$m18flag=0;$m19flag=0;$m20flag=0;$m21flag=0;$m22flag=0;$m23flag=0;
  $k++;
 }
}
