# If no data file option is given, describe script purpose and input.
if [ -z $1 ] || [ -z $2 ]
then
	echo;
	echo "# Command syntax is : 'sh FASTQ_trimming.sh [dataset_R1] [dataset_R2]'";
	echo "# ";
	echo "#        [dataset_R1] & [dataset_R2] : Left & right read files to be trimmed of unbalanced reads.";
	echo "# ";
	echo "# This file attempts to cleanup a pair of FASTQ files which are not validly formatted due to extra";
	echo "#      and/or unbalanced lines.  This formatting problem results in FASTQC crashing.";
	echo "# ";
	echo;
	exit 1;   # exit with error.
else
	findStr=".fastq";
	replaceStr=".residue.fastq";
	residueName1=$(echo $1 | sed -e "s/$findStr/$replaceStr/g");
	residueName2=$(echo $2 | sed -e "s/$findStr/$replaceStr/g");
	# If final files are found, don't process these data files.
	if [ -f $residueName1 ] || [ -f $residueName2 ]
	then
		echo "# ";
		echo "# Trimming of unbalanced read pairs has already been run on this datsaet.";
		echo "# ";
		exit 0;   # exit with no error.
	else
		## FASTQ format per line, repeating.
		# echo "@ id"
		# echo "sequence"
		# echo "+ id"
		# echo "quality"

		# number of lines of each file.
		echo -e "\tTotal number of lines in each file:";
		length_full1=$(wc -l $1 | awk '{print $1}');
		echo -e "\t\tFile1: "$length_full1;
		length_full2=$(wc -l $2 | awk '{print $1}');
		echo -e "\t\tFile2: "$length_full2;

		# modulus of the number of lines by 4, as the fastq format is in blocks of 4 lines.
		# this will give us the number of extra lines in each file.
		echo -e "\tLines at end of each file suggestive of a cropped entry:";
		length_extra1=$(expr $length_full1 % 4);
		echo -e "\t\tFile1: "$length_extra1;
		length_extra2=$(expr $length_full2 % 4);
		echo -e "\t\tFile2: "$length_extra2;

		# the number of properly formated lines per file.
		echo -e "\tNumber of properly formated lines per file:";
		length_base1=$(expr $length_full1 - $length_extra1);
		echo -e "\t\tFile1: "$length_base1;
		length_base2=$(expr $length_full2 - $length_extra2);
		echo -e "\t\tFile2: "$length_base2;

		# find the lesser length of properly formated lines per file.
		# this will be used to trim the original data files to a consistent valid length.
		echo -e "\tNumber of lines for paired-end reads in files:";
		stringer=$length_base1" "$length_base2;
		stringer=$(echo $stringer | tr ' ' '\n' | sort -n)
		lesser_length_base=$(echo $stringer | awk '{print $1}');
		echo -e "\t\t"$lesser_length_base;

		# find the residue lengths per file.
		# this is the number of lines of incomplete and extra reads.
		echo -e "\tNumber of lines at the end of each file to be removed:";
		length_residue1=$(expr $length_full1 - $lesser_length_base);
		echo -e "\t\tFile1: "$length_residue1;
		length_residue2=$(expr $length_full2 - $lesser_length_base);
		echo -e "\t\tFile2: "$length_residue2;

		findStr=".fastq";
		replaceStr=".trimmed.fastq";
		trimmedName1=$(echo $1 | sed -e "s/$findStr/$replaceStr/g");
		trimmedName2=$(echo $2 | sed -e "s/$findStr/$replaceStr/g");

		# make copies of the original read files which are trimmed to valid and the same length.
		echo -e "\tMaking trimmed valid fastq files:";
		echo -e "\t\tFile1: "$trimmedName1;
		cat $1 | head -$lesser_length_base > $trimmedName1;
		echo -e "\t\tFile2: "$trimmedName2;
		cat $2 | head -$lesser_length_base > $trimmedName2;

		# make residue files containing any incomplete or extra read lines.
		echo -e "\tMaking files containing residual lines:";
		if [ "$length_residue1" -gt 0 ]
		then
			echo -e "\t\tFile1: "$residueName1;
			cat $1 | tail -$length_residue1 > $residueName1;
		fi
		if [ "$length_residue2" -gt 0 ]
		then
			echo -e "\t\tFile2: "$residueName2;
			cat $2 | tail -$length_residue2 > $residueName2;
		fi

		echo -e "\tReplacing original files with trimmed versions.";
		# cleanup original raw files.
		rm $1;
		rm $2;
		# rename produced valid FASTQ files to original file names.
		mv $trimmedName1 $1
		mv $trimmedName2 $2
		echo;
	fi
fi
