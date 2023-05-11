input=$1
file="/home/mgray7/output2/fasta/$input.txt";
exec 5<$file
i=0;
while read line <&5; do

if (( $i%2==1 ))

then
# echo $line
echo $line > /home/mgray7/SparseRNAFolD/sparsemfefold/out1.txt
# /usr/bin/time -o out.txt -f "%e\t%M" ./build/SparseMFEFold -d2 $line;
/usr/bin/time -o out.txt -f "%e\t%M" /home/mgray7/includes/vienna/bin/RNAfold -d2 /home/mgray7/SparseRNAFolD/sparsemfefold/out1.txt;
# /usr/bin/time -o out.txt -f "%e\t%M" /home/mgray7/includes/sparsemfe/bin/SparseMFEFold $line;
# /usr/bin/time -o out.txt -f "%e\t%M" ../../LinearFold/linearfold -V -d 2 < /home/mgray7/SparseRNAFolD/sparsemfefold/out1.txt
cat "/home/mgray7/SparseRNAFolD/sparsemfefold/out.txt" >> "/home/mgray7/output2/results/$input.txt"
fi
#   
i=$((i+1));
done 
echo " first is $i";
exec 5</dev/null
# /usr/bin/time -o out.txt -f "%M\n%e" ./build/SparseMFEFold GGGAAACCC
# /usr/bin/time -o out.txt -f "%M\t%e" ./build/SparseMFEFold GGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCC