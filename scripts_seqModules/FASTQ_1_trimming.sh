# If no data file option is given, describe script purpose and input.
if [ -z $1 ]
then
	echo;
	echo "# Command syntax is : 'sh FASTQ_trimming.sh [dataset]'";
	echo "# ";
	echo "#        [dataset] : read file to be trimmed of incomplete final read.";
	echo "# ";
	echo "# This script attempts to cleanup a FASTQ file which is not validly formatted due to early";
	echo "#      termination of data transfer.  This formatting problem results in FASTQC crashing.";
	echo "# ";
	echo;
	exit 1;   # exit with error.
else
	findStr=".fastq";
	replaceStr=".residue.fastq";
	residueName1=$(echo $1 | sed -e "s/$findStr$/$replaceStr/g");
	# If final files are found, don't process these data files.
	if [ -f $residueName1 ]
	then
		echo;
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
		echo "#\tTotal number of lines in file:";
		length_full1=$(wc -l $1 | awk '{print $1}');
		echo "#\t\tFile1: "$length_full1;

		# modulus of the number of lines by 4, as the fastq format is in blocks of 4 lines.
		# this will give us the number of extra lines in each file.
		echo "#\tLines at end of each file suggestive of a cropped entry:";
		length_extra1=$(expr $length_full1 % 4);
		echo "#\t\tFile1: "$length_extra1;

		# If the number of extra lines is zero, then nothing needs to be done.
		# If the number of extra lines is more than zero, then create a trimmed file to replace the original.
		if [ "$length_extra1" -gt 0 ]
		then
			# the number of properly formated lines per file.
			echo "#\tNumber of properly formated lines per file:";
			length_base1=$(expr $length_full1 - $length_extra1);
			echo "#\t\tFile1: "$length_base1;

			# find the residue lengths per file.
			# this is the number of lines of incomplete and extra reads.
			echo "#\tNumber of lines at the end of file to be removed:";
			length_residue1=$length_extra1;
			echo "#\t\tFile1: "$length_residue1;

			findStr=".fastq";
			replaceStr=".trimmed.fastq";
			trimmedName1=$(echo $1 | sed -e "s/$findStr$/$replaceStr/g");

			# make copy of the original read file which are trimmed to valid length.
			echo "#\tMaking trimmed valid fastq files:";
			echo "#\t\tFile1: "$trimmedName1;
			cat $1 | head -$length_base1 > $trimmedName1;

			# make residue file containing any incomplete or extra read lines.
			echo "#\tMaking file containing residual lines:";
			if [ "$length_residue1" -gt 0 ]
			then
				echo "#\t\tFile1: "$residueName1;
				cat $1 | tail -$length_residue1 > $residueName1;
			fi

			echo "#\tReplacing original files with trimmed versions.";
			# cleanup original raw files.
			rm $1;
			# rename produced valid FASTQ files to original file names.
			mv $trimmedName1 $1
		else
			echo "#\tFile trimming does not need performed on this dataset.";
		fi
	fi
fi
