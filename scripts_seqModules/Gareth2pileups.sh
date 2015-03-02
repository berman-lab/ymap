#n/bash

#===================================================================================================================================
# Translate Gareth's pileup format into the pileup formats used by the pipeline.
#-----------------------------------------------------------------------------------------------------------------------------------
user=$1;
project=$2;
inputFile=$3;
main_dir=$(pwd)"/../../";

projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
tempDir=$projectDirectory"temp_dir/";
logFile=$projectDirectory"process_log.txt";

echo "projectDirectory   = '"$projectDirectory"'";
echo "tempDir            = '"$tempDir"'";
echo "logFile            = '"$logFile"'";


echo "\t\t| Shell : Convert the uploaded tab-delimited text data to pileup formats used in the pipeline." >> $logFile;
echo "\t\t|\t*===========================================================*" >> $logFile;
echo "\t\t|\t| Log of 'scripts_seqModules/Gareth2pileups.sh'             |" >> $logFile;
echo "\t\t|\t*-----------------------------------------------------------*" >> $logFile;
echo "\t\t|\t| Arguments:" >> $logFile;
echo "\t\t|\t|     user      = "$user >> $logFile;
echo "\t\t|\t|     project   = "$project >> $logFile;
echo "\t\t|\t|     main_dir  = "$main_dir >> $logFile;
echo "\t\t|\t|     inputFile = "$inputFile >> $logFile;
echo "\t\t|\t| " >> $logFile;
echo "\t\t|\t| inputFile is a tab-delimited-text file, with optional collumns in parentheses. " >> $logFile;
echo "\t\t|\t|         1\tChromosome name string." >> $logFile;
echo "\t\t|\t|         2\tBp coordinate." >> $logFile;
echo "\t\t|\t|         3\tPrimary base call." >> $logFile;
echo "\t\t|\t|         4\tNumber of reads with primary base call." >> $logFile;
echo "\t\t|\t|         (5)\tSecondary base call." >> $logFile;
echo "\t\t|\t|         (6)\tNumber of reads with secondary base call." >> $logFile;
echo "\t\t|\t|         (7)\tTertiary base call." >> $logFile;
echo "\t\t|\t|         (8)\tNumber of reads with tertiary base call." >> $logFile;
echo "\t\t|\t|         (9)\tQuaternary base call." >> $logFile;
echo "\t\t|\t|         (10)\tNumber of reads with quaternary base call." >> $logFile;
echo "\t\t|\t| " >> $logFile;
echo "\t\t|\t| Output files are placed in the project directory:" >> $logFile;
echo "\t\t|\t|        putative_SNPs_v4.txt" >> $logFile;
echo "\t\t|\t|        SNP_CNV_v1.txt" >> $logFile;
echo "\t\t|\t| " >> $logFile;
echo "\t\t|\t| Making temp dir in project folder." >> $logFile;

mkdir $tempDir;
outputFile1=$projectDirectory"putative_SNPs_v4.txt";
outputFile2=$projectDirectory"SNP_CNV_v1.txt";
tempFile1=$tempDir"temp.putative_SNPs_v4.txt";
tempFile2=$tempDir"temp.SNP_CNV_v1.txt";

echo "\t\t|\t| Processing txt file to produce pileups." >> $logFile;
# Loop through the file, line by line.
while read line
do
	set $line;

	# Rest values.
	chromosome="";
	coordinage="";
	allele1="";
	count1=0;
	allele2="";
	count2=0;
	allele3="";
	count3=0;
	allele4="";
	count4=0;

	# Read in new values from columns in line.
	chromosome=$1;
	coordinate=$2;
	allele1=$3;
	count1=$4;
	allele2=$5;
	count2=$6;
	allele3=$7;
	count3=$8;
	allele4=$9;
	count4=$10;

	A_count=0;
	T_count=0;
	G_count=0;
	C_count=0;
	existA=0;
	existT=0;
	existG=0;
	existC=0;

	if [ "$allele1" = "A" ]; then
		A_count=$count1;
		existA=1;
	elif [ "$allele1" = "T" ]; then
		T_count=$count1;
		existT=1;
	elif [ "$allele1" = "G" ]; then
		G_count=$count1;
		existG=1;
	elif [ "$allele1" = "C" ]; then
		C_count=$count1;
		existC=1;
	fi

	if [ "$allele2" = "A" ]; then
		A_count=$count2;
		existA=1;
	elif [ "$allele2" = "T" ]; then
		T_count=$count2;
		existT=1;
	elif [ "$allele2" = "G" ]; then
		G_count=$count2;
		existG=1;
	elif [ "$allele2" = "C" ]; then
		C_count=$count2;
		existC=1;
	fi

	if [ "$allele3" = "A" ]; then
		A_count=$count3;
		existA=1;
	elif [ "$allele3" = "T" ]; then
		T_count=$count3;
		existT=1;
	elif [ "$allele3" = "G" ]; then
		G_count=$count3;
		existG=1;
	elif [ "$allele3" = "C" ]; then
		C_count=$count3;
		existC=1;
	fi

	if [ "$allele4" = "A" ]; then
		A_count=$count4;
		existA=1;
	elif [ "$allele4" = "T" ]; then
		T_count=$count4;
		existT=1;
	elif [ "$allele4" = "G" ]; then
		G_count=$count4;
		existG=1;
	elif [ "$allele4" = "C" ]; then
		C_count=$count4;
		existC=1;
	fi
	total_count=$(($A_count + $T_count + $G_count + $C_count));

	# Line from pipeline processing script defining 'putative_SNPs_v4.txt' file format.
	#	print chrom + '\t' + pos + '\t' + ref_base + '\t' + str(A) + '\t' + str(T) + '\t' +  str(G) + '\t' +  str(C)
	#
	# Line from pipeline processing script defining 'SNP_CNV_v1.txt' file format.
	#	print chrom + '\t' + pos + '\t' + str(total) + '\t' + ref_base + '\t' + str(A) + '\t' + str(T) + '\t' +  str(G) + '\t' +  str(C);
	#

	# Output line to 'putative_SNPs_v4.txt' file.
	if [ $(($existA+$existT+$existG+$existC)) -gt 1 ]; then
		output=$chromosome"\t"$coordinate"\t"$allele1"\t"$A_count"\t"$T_count"\t"$G_count"\t"$C_count;
		echo $output >> $tempFile1;
	fi

	# Output line to 'SNP_CNV_v1.txt' file.
	output=$chromosome"\t"$coordinate"\t"$total_count"\t"$allele1"\t"$A_count"\t"$T_count"\t"$G_count"\t"$C_count;
	echo $output >> $tempFile2;


	# Let the user know the script is proceding through chromosomes...
	if [ "$previous_chr" != "$chromosome" ]; then
		echo "Current chromosome = "$chromosome;
	fi
	previous_chr=$chromosome;
done < $projectDirectory$inputFile

echo "\t\t|\t| Moving temp files to final location." >> $logFile;
mv $tempFile1 $outputFile1;
mv $tempFile2 $outputFile2;

rm -rf $tempDir;

echo "\t\t|\t*-----------------------------------------------------------*" >> $logFile;
echo "\t\t|\t| 'scripts_seqModules/Gareth2pileups.sh' has completed.     |" >> $logFile;
echo "\t\t|\t*===========================================================*" >> $logFile;
