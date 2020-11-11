open IN1, "<$ARGV[0]";

$mut_block=0;
$json_depth=0;
$mut_block2=0;
while(<IN1>)
{
        $line=$_;
        #print "$json_depth | $mut_block | $line";
        # #if($mut_clock)

        # #if(/\"nodes\"\:/){$mut_block=1;}
        # #if( $mut_block == 1)
        # {
        #       if(/\"(.+)\"/ && $line !~ /\"sequence\"\:/)     #warning: you should not have sequence id "sequence"
        #       {
        #               #print "$node_id\t$mut_cnt\n";
        #               $node_id=$1; 
        #               $mut_cnt=0;
        #               #print "$node_id\n";
        #               while(<IN1>)
        #               {
        #                       if(/\]/){last;}
        #                       $mut_cnt++;
        #               }
        #               print "!$node_id\t$mut_cnt\n";
        #       }

        # }

        if( $mut_block == 1 && $json_depth==2)
        {
                /\"(.+)\"/; 
                $json_depth+=1;
                #print "!$1\n";
                $node_id=$1; 
                $mut_cnt=0;
                $del_cnt=0;
                while(<IN1>)
                {
                        $line=$_;
                        #print "?$_";
                        if($line=~/]/) {last;}
                        if($line!~/-/){$mut_cnt++;}
                        else{$del_cnt=1;}
                }
                if($mut_cnt > 0){$mut_cnt-=1;}
                $mut_cnt+=$del_cnt;
                $mut_cnt_hash{$node_id}=$mut_cnt/29903;
        }

        if(/\{/){$json_depth+=1;}
        if(/\}/){$json_depth-=1;}
        if(/\"nodes\"\:/){$mut_block=1;}
        if($mut_block == 1 && $json_depth == 1 ){$mut_block=0;}
        #print "$json_depth | $mut_block | $line";
}
#print "\n";

open IN2, "<$ARGV[1]";

while(<IN2>)
{
        chomp;
        @array=split /:/;
        $str2=$array[0];

        for($loop=1;$loop<=$#array;$loop++)
        {
                $array[$loop-1]=~/([^\(\)\,]+)$/;
                $read_id=$1;
                $array[$loop]=~/([\(\)\,].+)/;
                $length=$1;
                $str2=$str2.":".$mut_cnt_hash{$read_id}.$length;
                #$str=$str2.$mut_cnt_hash{$read_id};
                #print "$str\n";
        }

        foreach ( @array )
    {
        $element=$_;
                /([^\(\)\,]+)$/;
                $read_id=$1;
                /([\(\)\,].+)/;
                $length=$1;
                #print "$element\t$length\t$mut_cnt_hash{$read_id}\n";
        }
        #print "eol\n";
        print "$str2\n";
}
