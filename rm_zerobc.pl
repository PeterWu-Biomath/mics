open IN1,"gzip -dc $ARGV[0] | ";
open OUT,"| gzip > $ARGV[1]";

while(<IN1>)
{
        $S=$_;
        if(/$S=~0_0_0/)
        {
                <IN1>;
                <IN1>;
                <IN1>;
        }
        else
        {
                print OUT "$S";
                $S=<IN1>;
                print OUT "$S";
                $S=<IN1>;
                print OUT "$S";
                $S=<IN1>;
                print OUT "$S";
        }
}
