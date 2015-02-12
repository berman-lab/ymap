##---------------------------------------------------------------------------------------------
## Reformat fasta : puts FASTA sequences into single lines, then sorts entries by seq length.
##---------------------------------------------------------------------------------------------
echo "[[=- Reformatting FASTA entries into single lines per sequence. -=]]";
# 0) Called like : "sh reformat_fasta.sh file.fa"
# 1) For lines that start with ">", convert the ending "\n" into "\t".
# 2) Removes all newline characters.
# 3) Adds newlines in front of ">"s to split FASTA entries onto separate lines.
# 4) Replaces the tab characters with newlines, to restore proper fasta format.
# 5) Remove initial blank lines...  added above as artifact of adding lines between entries.

# 1) For lines that start with ">", convert the ending "\n" into "\t".
sed -i '/^[>]/{N;s/\n/\t/}' $1

# 2) Removes all newline characters.
perl -pi -e 's/\n//g' $1

# 3) Adds newlines in front of ">"s to split FASTA entries onto separate lines.
perl -pi -e 's/>/\n>/g' $1

# 4) Replaces the tab characters with newlines, to restore proper fasta format.
perl -pi -e 's/\t/\n/g' $1

# 5) Remove initial blank lines...  added above as artifact of adding lines between entries.
awk 'NR > 1 { print }' < $1 > $1.temp
rm $1;
mv $1.temp $1;

##---------------------------------------------------------------------------------------------
## Reformat fasta : computes reverse complement of each FASTA entry.
##---------------------------------------------------------------------------------------------
echo "[[=- Determining reverse-complement for sequence in each FASTA entry. -=]]";
outFile=$1.temp;
count=0;
while read p; do
        count=$(expr $count + 1);
        headerTest=${p:0:1};
        if [ "$headerTest" == ">" ]
        then
                echo $p >> $outFile;
                echo -e "\theader line.";
        else
                lineLength=${#p}
                if [ $lineLength > 1 ]
                then
                        echo -e "\tSequence line";
                else
                        echo -e "\tBlank line";
                fi
                # Generate reverse-complement versions of the sequence input.
                seq=$p;
                seq_length_1=${#seq};
                seq1RC="";
                for ((i=$seq_length_1-1; i>=0; --i))
                do
                        nt1=${seq:$i:1};
                        # Find the complement of the base in compared strand 1.
                        if [ "$nt1" == "A" ]; then nt1_="T"; fi;        if [ "$nt1" == "a" ]; then nt1_="t"; fi;
                        if [ "$nt1" == "T" ]; then nt1_="A"; fi;        if [ "$nt1" == "t" ]; then nt1_="a"; fi;
                        if [ "$nt1" == "G" ]; then nt1_="C"; fi;        if [ "$nt1" == "g" ]; then nt1_="c"; fi;
                        if [ "$nt1" == "C" ]; then nt1_="G"; fi;        if [ "$nt1" == "c" ]; then nt1_="g"; fi;
                        if [ "$nt1" == "Y" ]; then nt1_="R"; fi;        if [ "$nt1" == "y" ]; then nt1_="r"; fi;
                        if [ "$nt1" == "R" ]; then nt1_="Y"; fi;        if [ "$nt1" == "r" ]; then nt1_="y"; fi;
                        if [ "$nt1" == "S" ]; then nt1_="S"; fi;        if [ "$nt1" == "s" ]; then nt1_="s"; fi;
                        if [ "$nt1" == "W" ]; then nt1_="W"; fi;        if [ "$nt1" == "w" ]; then nt1_="w"; fi;
                        if [ "$nt1" == "K" ]; then nt1_="M"; fi;        if [ "$nt1" == "k" ]; then nt1_="m"; fi;
                        if [ "$nt1" == "M" ]; then nt1_="K"; fi;        if [ "$nt1" == "m" ]; then nt1_="k"; fi;
                        if [ "$nt1" == "B" ]; then nt1_="V"; fi;        if [ "$nt1" == "b" ]; then nt1_="v"; fi;
                        if [ "$nt1" == "V" ]; then nt1_="B"; fi;        if [ "$nt1" == "v" ]; then nt1_="b"; fi;
                        if [ "$nt1" == "D" ]; then nt1_="H"; fi;        if [ "$nt1" == "d" ]; then nt1_="h"; fi;
                        if [ "$nt1" == "H" ]; then nt1_="D"; fi;        if [ "$nt1" == "h" ]; then nt1_="d"; fi;
                        if [ "$nt1" == "N" ]; then nt1_="N"; fi;        if [ "$nt1" == "n" ]; then nt1_="n"; fi;
                        nt1=$nt1_;  
                        seq1RC=$seq1RC$nt1_;
                done
                echo $seq1RC >> $outFile;
        fi
done < $1
rm $1;
mv $1.temp $1;

##---------------------------------------------------------------------------------------------
## Reformat fasta : inserts newline characters into long sequences to break into 100bp lines.
##---------------------------------------------------------------------------------------------
echo "[[=- Reformatting FASTA entries to wrap sequence at 100 bp. -=]]";
# 0) Called like : "sh reformat_fasta.sh file.fa"
# 1) Adds a newline between each fasta entry, for ease of reading.
# 2) Adds a newline between each fasta entry, for ease of reading.
# 3) Remove initial blank lines...  added above as artifact of adding lines between entries.

# 1) Adds a newline between each fasta entry, for ease of reading.
perl -pi -e 's/>/\n>/g' $1

# 2) Adds a newline between each fasta entry, for ease of reading.
perl -pi -e 's/.{100}/$&\n/g' $1

# 3) Remove initial blank lines...  added above as artifact of adding lines between entries.
awk 'NR > 1 { print }' < $1 > $1.temp
rm $1;
mv $1.temp $1;
