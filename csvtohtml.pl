open IN1,"</home/jovyan/scripts/csvtemplate.html";

while(<IN1>){print;}

close IN1;
open IN1,"<$ARGV[0]";

$line=<IN1>;

@array=split /\t/ ,$line;

print "<thead><tr>";
foreach $e (@array){print "<th>$e</th>";}
print "</tr></thead>\n";

while($line= <IN1>)
{
	@array=split /\t/ ,$line;
	print "<tbody><tr>";
	foreach $e (@array){print "<td>$e</td>";}
	print "</tr></tbody>\n";
}

print "</table></div>	</div></div></body></html>\n";

close IN1;
