##---------------------------------------------------------------------------------------------
## Reformat fasta 2 : converts single-line FASTA entries into multi-line entries.
##---------------------------------------------------------------------------------------------
# 0) Called like : "sh FASTA_reformat_2.sh file.fa"
# 1) Adds newlines in front of ">"s to split FASTA entries onto separate lines.
# 2) Adds a newline after every 100bp of sequence.
# 3) Remove initial blank lines...  added above as artifact of adding lines between entries.


# If no data file option is given, describe script purpose and input.
if [ -z $1 ] 
then
	echo;
	echo "# Command syntax is : 'sh FASTQ_reformtat_1.sh [FASTA seq file]'";
	echo "# ";
	echo "#        [FASTA seq file] : Genome sequence file in FASTA format.";
	echo "# ";
	echo "# This script will take a file containing single-line FASTA entries and reformat";
	echo "# them to have one header line and many sequence lines 100bp long per entry.";
	echo;
	exit 1;
else
	# 1) Adds newlines in front of ">"s to split FASTA entries onto separate lines.
	perl -pi -e 's/>/\n>/g' $1

	# 1) For lines that don't start with ">", add a newline after every 100 characters.
	perl -pi -e 'if (!/^[>]/) { s/.{100}/$&\n/g }' $1;

	# 3) Remove initial blank line added above as artifact of adding lines between entries.
	awk 'NR > 1 { print }' < $1 > $1.temp
	rm $1;
	mv $1.temp $1;
fi
