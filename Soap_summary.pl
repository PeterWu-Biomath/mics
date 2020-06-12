open IN1,"<$ARGV[0]/result/02.QCL/test/Basic_Statistics_of_Sequencing_Quality.txt";
$line_cnt=0;
while(<IN1>)
{
        chomp;
        $line=$_;
        @array=split /\s+/ ,$line;
        #print "$line_cnt\t";
        #$cnt=0;
        #foreach $x (@array){print "$cnt:$x\t";$cnt++;}
        #print "\n";
        if($line_cnt==1){$read_length=$array[2];}
        if($line_cnt==3){$filtered_read=$array[5]*2;}
        if($line_cnt==2){$raw_read_cnt=$array[4]*2;}
        if($line_cnt==16)
        {
                $raw_read_Q20=$array[13]+$array[17];
                $clean_read_Q20=$array[15]+$array[19];
        }
        if($line_cnt==17)
        {
                $raw_read_Q30=$array[13]+$array[17];
                $clean_read_Q30=$array[15]+$array[19];
        }
        if($line_cnt==10)
        { 
                $raw_read_A=$array[5]+$array[9];
                $clean_read_A=$array[7]+$array[11];
        }
        if($line_cnt==11)
        { 
                $raw_read_C=$array[5]+$array[9];
                $clean_read_C=$array[7]+$array[11];
        }
        if($line_cnt==12)
        { 
                $raw_read_G=$array[5]+$array[9];
                $clean_read_G=$array[7]+$array[11];
        }
        if($line_cnt==13)
        { 
                $raw_read_T=$array[5]+$array[9];
                $clean_read_T=$array[7]+$array[11];
        }
        $line_cnt++;
}

close IN1;

$raw_base=$raw_read_cnt*$read_length;
$clean_reads=$raw_read_cnt-$filtered_read;
open IN1,"<$ARGV[0]/result/02.QCL/test/Statistics_of_Filtered_Reads.txt";
$line_cnt=0;
while(<IN1>)
{
        chomp;
        $line=$_;
        @array=split /\s+/ ,$line;
        #print "$line_cnt\t";
        #$cnt=0;
        #foreach $x (@array){print "$cnt:$x\t";$cnt++;}
        #print "\n";
        if($line_cnt==2){$adaptor_cnt=$array[4];}
        if($line_cnt==5){$dup_cnt=$array[4];}
        if($line_cnt==6){$N_cnt=$array[6];}
        if($line_cnt==3){$low_qual_cnt=$array[5];}
        $line_cnt++;
}

$adaptor_ratio=int($adaptor_cnt*10000/$raw_read_cnt+0.5)/100;
$N_ratio=int($N_cnt*10000/$raw_read_cnt+0.5)/100;
$low_qual_ratio=int($low_qual_cnt*10000/$raw_read_cnt+0.5)/100;
$dup_ratio=int($dup_cnt*10000/$raw_read_cnt+0.5)/100;
$clean_base=$clean_reads*$read_length;
$raw_read_Q20_ratio=int($raw_read_Q20*10000/$raw_base+0.5)/100;
$raw_read_Q30_ratio=int($raw_read_Q30*10000/$raw_base+0.5)/100;
$clean_read_Q20_ratio=int($clean_read_Q20*10000/$clean_base+0.5)/100;
$clean_read_Q30_ratio=int($clean_read_Q30*10000/$clean_base+0.5)/100;
$clean_GC=int(10000*($clean_read_G+$clean_read_C)/($clean_read_T+$clean_read_A+$clean_read_G+$clean_read_C)+0.5)/100;
$raw_GC=int(10000*($raw_read_G+$raw_read_C)/($raw_read_T+$raw_read_A+$raw_read_G+$raw_read_C)+0.5)/100;

close IN1;
print "$ARGV[1]\t$read_length\t$raw_read_cnt\t$raw_base\t";
print "$N_cnt\t$N_ratio%\t$low_qual_cnt\t$low_qual_ratio%\t$dup_cnt\t$dup_ratio%\t$adaptor_cnt\t$adaptor_ratio%\t";
print "$clean_reads\t$clean_bases\t$raw_read_Q20_ratio%\t$clean_read_Q20_ratio%\t$raw_read_Q30_ratio%\t$clean_read_Q30_ratio%\t$raw_GC%\t$clean_GC%\n";
