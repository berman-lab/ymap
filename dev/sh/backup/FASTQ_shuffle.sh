# If no data file option is given, describe script purpose and input.
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]
then
	echo;
	echo "# Command syntax is : 'sh FASTQ_shuffle.sh [dataset_R1] [dataset_R2] [outfile]'";
	echo "# ";
	echo "#        [dataset_R1] & [dataset_R2] : Left & right read files to be shuffled into one FASTQ file.";
	echo "#        [outfile]                   : Destination file for FASTQ entries.";
	echo "# ";
	echo "# This script expects the two input files to have the same number of lines.";
	echo "# ";
	echo;
	exit 1;
else
	## FASTQ format per line, repeating.
	# echo "@ id";
	# echo "sequence"
	# echo "+ id"
	# echo "quality"

	# number of lines of each file.
	length_full1=$(wc -l $1 | awk '{print $1}');
	length_full2=$(wc -l $2 | awk '{print $1}');

	if [ $length_full1 -ne $length_full2 ]
	then
		echo "##";
		echo "## Input FASTQ files have different numbers of lines.";
		echo "##";
		exit;
	else
		# initialize and clear output file.
		echo "null" >> $3;
		cp /dev/null $3;

		# open extra file descriptors for input.
		exec 4< $1;
		exec 5< $2;

		# read lines from input files, then append to output file in blocks of four lines.
		while read line1_1 <&4;
		do
			read line1_2 <&4;
			read line1_3 <&4;
			read line1_4 <&4;
			read line2_1 <&5;
			read line2_2 <&5;
			read line2_3 <&5;
			read line2_4 <&5;

			echo $line1_1 >> $3;
			echo $line1_2 >> $3;
			echo $line1_3 >> $3;
			echo $line1_4 >> $3;
			echo $line2_1 >> $3;
			echo $line2_2 >> $3;
			echo $line2_3 >> $3;
			echo $line2_4 >> $3;
		done

		# close the extra file descriptors.
		exec 4<&-;
		exec 5<&-;
	fi
fi
