##---------------------------------------------------------------------------------------------
## Reformat fasta 1 : converts multi-line FASTA entries into single-line entries.
##---------------------------------------------------------------------------------------------
# 0) Called like : "sh FASTA_reformat_1.sh file.fa"
# 1) For lines that start with ">", convert the ending "\n" into "\t".
# 2) Removes all newline characters.
# 3) Adds newlines in front of ">"s to split FASTA entries onto separate lines.
# 4) Replaces the tab characters with newlines, to restore proper fasta format.
# 5) Remove initial blank lines...  added above as artifact of adding lines between entries.


# If no data file option is given, describe script purpose and input.
if [ -z $1 ]
then
	echo;
	echo "# Command syntax is : 'sh FASTA_reformat_1.sh [FASTA seq file]'";
	echo "# ";
	echo "#        [FASTA seq file] : Genome sequence file in FASTA format.";
	echo "# ";
	echo "# This script will take a file containing multi-line FASTA entries and reformat";
	echo "# them to have only one line for the header and for the sequence for each entry."; 
	echo "# ";
	echo;
	exit 1;
else
	# 1) For lines that start with ">", convert the ending "\n" into "\t".
	cp $1 $1.temp.bak
	TAB=$'\t'                                                           # Mac's sed is really old and doesn't understan "\t", thus we must resort to trickery. Also not giving -i an extesion doesn't work for some reason.
	sed -e '/^[>]/{' -e "N;s/\n/${TAB}/" -e '}' -i $1.temp.bak $1
	rm $1.temp.bak

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
fi
