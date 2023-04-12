input=$1
file="/home/mgray7/output2/fasta/$input.txt";
file2="/home/mgray7/output2/structures/$input.txt";
exec 5<$file
exec 6<$file2
while read line1 <&5 && read line2 <&6; do

if (( $i%2==1 ))

then
/usr/bin/time -o out.txt -f "%e\t%M" ./build/SparseMFEFoldTriplet -d1 $line1;
cat "/home/mgray7/sparsemfe/sparsemfefold/out.txt" >> "/home/mgray7/output2/results/$input.txt"
fi
i=$((i+1));
done
echo "first is $i";
#   echo "${line}";
exec 5</dev/null
exec 6</dev/null


# echo $line > /home/mgray7/sparsemfe/sparsemfefold/out1.txt
# /usr/bin/time -o out.txt -f "%e\t%M" ../../LinearFold/linearfold -V -d 2 < /home/mgray7/sparsemfe/sparsemfefold/out1.txt

#   

# done < "${file}"

