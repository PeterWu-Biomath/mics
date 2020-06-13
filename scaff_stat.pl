my $file1 =  $ARGV[0];

open IN1, "<$file1";

while(<IN1>)
{
        chomp;
        if(/>/){push @scaffs,"";}
        else{$scaffs[-1]=$scaffs[-1].$_;}
}

my @scaff_len;
my @contig_len;
my $GC_cnt=0;
my $N_cnt=0;
my $AT_cnt=0;

foreach my $scaff (@scaffs) {
            push @scaff_len, length($scaff);
            foreach $char (split //, $scaff) {
                      if($char eq "N" || $char eq "n"){$N_cnt++;}
                      elsif($char eq "A" || $char eq "a" || $char eq "T" || $char eq "t"){$AT_cnt++;}
                      elsif($char eq "C" || $char eq "c" || $char eq "G" || $char eq "g"){$GC_cnt++;}
              }
              foreach $contig (split /n|N/,$scaff )
              {
                if(length($contig)>0){push @contig_len, length($contig);}
              }
        }
    
my @sort_scaff=sort{$b <=> $a}(@scaff_len);
my @sort_contig=sort{$b <=> $a}(@contig_len);
#foreach my $contig_lenx (@sort_contig) {print "$contig_lenx\n";}

$sample= $ARGV[1];

$genome_size=$GC_cnt+$N_cnt+$AT_cnt;

$GC_ratio=int(($GC_cnt*10000)/ ($GC_cnt+$AT_cnt) +0.5)/100;

$scaffold_num=@scaff_len;

$acc=0;
foreach my $scaff_len (@sort_scaff)
{
        $acc=$acc+$scaff_len;
        if($acc*2 > $genome_size){$scaff_N50=$scaff_len; last;}
}

$acc=0;
foreach my $scaff_len (@sort_scaff)
{
        $acc=$acc+$scaff_len;
        if($acc*10 > $genome_size*9){$scaff_N90=$scaff_len; last;}
}

$average_scaff=int(100*$genome_size/$scaffold_num)/100;
$max_scaff=$sort_scaff[-1];
$min_scaff=$sort_scaff[0];

$contig_num=@sort_contig;

$acc=0;
$contig_N50=0;
foreach my $scaff_len (@sort_contig)
{
        $acc=$acc+$scaff_len;
        #       print $acc*2;
        #print "\t";
        #print ($GC_cnt+$AT_cnt);
        #print "\n";
        if($acc*2 > ($GC_cnt+$AT_cnt)){$contig_N50=$scaff_len; last;}
}
#print "$contig_N50\n";
$acc=0;
foreach my $scaff_len (@sort_contig)
{
        $acc=$acc+$scaff_len;
        if($acc*10 > ($GC_cnt+$AT_cnt)*9){$contig_N90=$scaff_len; last;}
} 

$average_contig=int(100*($GC_cnt+$AT_cnt)/$contig_num)/100;

$max_contig=$sort_contig[-1];
$min_contig=$sort_contig[0];


print "$sample\t$genome_size\t$GC_ratio%\t$scaffold_num\t$scaff_N50\t$scaff_N90\t$average_scaff\t$max_scaff\t$min_scaff\t";
print "$contig_num\t$contig_N50\t$contig_N90\t$average_contig\t$max_contig\t$min_contig\t$N_cnt\n"
