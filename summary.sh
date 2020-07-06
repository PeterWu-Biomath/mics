echo "Sample	Total Reads Unclassified Reads Number	Unclassified Rate	Classified Reads Number	Classified Rate	Genome Size" > kraken.csv
echo "Sample	Read Length	Raw_Reads	Raw_Base	Reads_with_N	Reads_with_N_Rate	Reads_with_LowQuality	Reads_with_LowQuality_Rate	Reads_with_Adapter	Reads_with_Adapter_Rate	Reads_with_Duplications	Reads_with_Duplications_Rate	SOAPnuke_Clean_Reads	SOAPnuke_Clean_Base	SOAPnuke_Clean_Base_Rate	Raw_Base_Q20	Clean_Base_Q20	Raw_Base_Q30	Clean_Base_Q30	GC_of_Raw_Base	GC_of_Clean_Base" > Soap.csv
echo "Sample	Genome size	GC content	Scaffold num	Scaffold N50 (bp)	Scaffold N90 (bp)	Scaffold average length (bp)	Scaffold maximum length (bp)	Scaffold minimum length (bp)	Contig num	Contig N50 (bp)	Contig N90 (bp)	Contig average length (bp)	Contig maximum length (bp)	Contig minimum length (bp)	Gap number (bp)
" > asm.csv

for dir in ../output/COV*/
do 
	dir1=${dir%*/}
	dir2=${dir1##*/}
       	perl kraken_summary.pl $dir $dir2 >> kraken.csv
       	perl Soap_summary.pl $dir $dir2 >> Soap.csv
       	perl asm_summary.pl $dir $dir2 >> asm.csv
done

perl csvtohtml.pl kraken.csv > ../output/kraken.html
perl csvtohtml.pl Soap.csv > ../output/soap.html
perl csvtohtml.pl asm.csv > ../output/asm.html
