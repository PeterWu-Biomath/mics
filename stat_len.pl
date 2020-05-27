$len=0;
while(<>)
{
        if(/\>/)
        {
                if($len > 0){print "$len\n";}
                $len=0;
        }
        else{$len=$len+length($_)-1;}
}
