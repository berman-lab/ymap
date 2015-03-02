# If no data file option is given, describe script purpose and input.
if [ -z $1 ] || [ -z $2 ]
then
	echo;
	echo "# Command syntax is : 'sh FASTQ_chop.sh [dataset] [length]'";
	echo "# ";
	echo "#        [dataset] : File containing FASTQ entries.";
	echo "#        [length]  : Number of lines at which to cut file in two.";
	echo "# ";
	echo "# This function is useful in trouble-shooting FASTQ file formatting issues.";
	echo "# ";
	echo;
	exit 1;
else
	length_full=$(wc -l $1 | awk '{print $1}');
	end_length=$(( $length_full - $2 ));
	head $1 -n $2 > ${1//.fastq/.1.fastq};
	tail $1 -n $end_length > ${1//.fastq/.2.fastq};
	rm $1;
fi
