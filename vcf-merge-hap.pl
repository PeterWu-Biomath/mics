#!/usr/bin/perl
# Author: Wu Pei

use strict;
use warnings;


my $vcf_file1 =  $ARGV[0];
my $vcf_file2 =  $ARGV[1];
open IN1, "<$vcf_file1";
open IN2, "<$vcf_file2";

print "#CHROM\tPOS\tREF\tALT\tHAP\n";

my (@chr1, @pos1,@ref1,@alt1);
my (@chr2, @pos2,@ref2,@alt2);
my $chr_num=0;
my %chr_hash;
while(<IN1>)
{
    chomp;  # avoid \n on last field
    my @array = split;
   next if($array[0] eq "R");
   if(exists ($chr_hash{$array[1]}))
   {
    push @chr1, $chr_hash{$array[1]};
   }
  else
  {
  $chr_num++;
  $chr_hash{$array[1]}=$chr_num;
  
    push @chr1, $chr_hash{$array[1]};
  }
   push @pos1, ($array[2]+0);
   push @ref1, $array[6];
   push @alt1, $array[7];
}

while(<IN2>)
{
    chomp;  # avoid \n on last field
    my @array = split;
   next if($array[0] eq "R") ;
   if(exists ($chr_hash{$array[1]}))
   {
    push @chr2, $chr_hash{$array[1]};
   }
  else
  {
  $chr_num++;
  $chr_hash{$array[1]}=$chr_num;
  
    push @chr2, $chr_hash{$array[1]};
  }
   push @pos2, ($array[2]+0);
   push @ref2, $array[6];
   push @alt2, $array[7];
}
my @chr_list;
while (my ($key, $value) = each (%chr_hash))
{
  $chr_list[$value] = $key;
}
my ($i,$j);
$i=0;
$j=0;

while( $i <= $#chr1 && $j<=$#chr2)
{
    my $cmp1=($chr1[$i] <=> $chr2[$j] );
    my $cmp2=($pos1[$i] <=> $pos2[$j] );
    my $cmp3=($ref1[$i] cmp $ref2[$j] );
    my $cmp4=($alt1[$i] cmp $alt2[$j] );
    my $cmp= ($cmp1*9+$cmp2*3+$cmp3);
    if($cmp>0)
    {
        print "$chr_list[$chr2[$j]]\t$pos2[$j]\t$ref2[$j]\t$alt2[$j]\t0|1\n";
        $j++;
    }
    elsif($cmp<0)
    {
        print "$chr_list[$chr1[$i]]\t$pos1[$i]\t$ref1[$i]\t$alt1[$i]\t1|0\n";
    $i++;
    }
    else
    {
    if($cmp4==0)
    {
        
        print "$chr_list[$chr1[$i]]\t$pos1[$i]\t$ref1[$i]\t$alt1[$i]\t1|1\n";
    }
    else
    {
    print "$chr_list[$chr1[$i]]\t$pos1[$i]\t$ref1[$i]\t$alt1[$i],$alt2[$j]\t1|2\n";
    }
    $i++;$j++;
    }
}
close IN1;
close IN2;