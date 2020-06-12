$sample =$ARGV[1];

open IN1,"<$ARGV[0]/result/01.KRAKEN/all_classify_kraken.xls";

<IN1>;
$line=<IN1>;
 close IN1;
 @array = split /\W+/, $line;
 #print $line;
 $total_reads=$array[1];
 $cls_reads=$array[5];
 $uncls_reads= $array[1]-$array[5];
 $uncls_ratio= int(10000*$uncls_reads/$total_reads)/100;
 $cls_ratio= int(10000*$cls_reads/$total_reads)/100;
open IN2, "<$ARGV[0]/result/05.ASS/test/scaffolds.fasta.stat";

<IN2>;
<IN2>;
<IN2>;
<IN2>;
$line=<IN2>;
close IN2;
@array = split /\W+/, $line;
$genome_size=$array[4];


 print "$sample\t$total_reads\t$uncls_reads\t$uncls_ratio%\t$cls_reads\t$cls_ratio%\t$genome_size\n";
