##---------------------------------------------------------------------------------------------
## Count subseq : counts the number of times a sequence string is found in a file.
##---------------------------------------------------------------------------------------------
# Called like : "sh FASTA_count_subseq.sh [string] [file.fa]"

# If no data file option is given, describe script purpose and input.
if [ -z $1 ]
then
	echo;
	echo "# Command syntax is : 'sh FASTA_count_subseq.sh [sequence] [FASTA seq file]'";
	echo "# ";
	echo "#        [sequence]       : Short DNA sequence string.";
	echo "#        [FASTA seq file] : Genome sequence file in FASTA format.";
	echo "# ";
	echo "# This script will count the incidence of the target sequence in the FASTA file.";
	echo "# The FASTA file should be reformated to single-line entries using FASTA_reformat_1.sh first.";
	echo "# ";
	echo;
	exit 1;
else
	grep -o $1 $2 | wc -l;
fi
