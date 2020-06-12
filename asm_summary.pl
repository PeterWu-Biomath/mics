$sample=$ARGV[1];

open IN1,"<$ARGV[0]/result/05.ASS/test/scaffolds.fasta.stat";
$line_cnt=0;
while(<IN1>)
{
        chomp;
        $line=$_;
        $line=~/:(.+)/;
        @array={};
        @array=split /\s+/ ,$1;
        if($line_cnt==3) {$scaffold_num=$array[1];$contig_num=$array[2];}
        if($line_cnt==4) {$genome_size=$array[1];}
        if($line_cnt==5) {$N_cnt=$array[1];}
        if($line_cnt==6) {$average_scaff=$array[1];$average_contig=$array[2];}
        if($line_cnt==7) {$scaff_N50=$array[1];$contig_N50=$array[2];}
        if($line_cnt==8) {$scaff_N90=$array[1];$contig_N90=$array[2];}
        if($line_cnt==9) {$max_scaff=$array[1];$max_contig=$array[2];}
        if($line_cnt==10) {$min_scaff=$array[1];$min_contig=$array[2];}
        if($line_cnt==11) {$GC_ratio=$array[2];}
        #print "$line_cnt\t$array[1]\t$array[2]\n";
        $line_cnt++;
}



print "$sample\t$genome_size\t$GC_ratio\t$scaffold_num\t$scaff_N50\t$scaff_N90\t$average_scaff\t$max_scaff\t$min_scaff\t";
print "$contig_num\t$contig_N50\t$contig_N90\t$average_contig\t$max_contig\t$min_contig\t$N_cnt\n";
