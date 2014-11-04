#===================================================================================================================================
# Decompose Sam files into FASTQ files for introduction into the sequence analysis pipeline.
#-----------------------------------------------------------------------------------------------------------------------------------
user=$1;
project=$2;
main_dir=$3;
inputFile=$4;
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logFile=$projectDirectory"process_log.txt";
echo "#|---- sam2fastw.sh ---- begin." >> $logFile;

# import locations of auxillary software for pipeline analysis.
. $main_dir"sh/local_installed_programs.sh";

# #===================================================================================================================================
# # Run picard-tools SamToFastq : This method requires large amounts of memory and will crash when it runs out of memory.
# #-----------------------------------------------------------------------------------------------------------------------------------
# java -Xmx16g -jar $picardDirectory"SamToFastq.jar" INPUT=$projectDirectory$inputFile FASTQ=$projectDirectory"data_r1.a.fastq" SECOND_END_FASTQ=$projectDirectory"data_r2.a.fastq";


#===================================================================================================================================
# Alternate method using generalized commandline tools with low memory footprint.
#-----------------------------------------------------------------------------------------------------------------------------------
echo "#| Translating SAM/BAM file to FASTQ format." >> $logFile;
echo "#| input file name             = "$inputFile >> $logFile;
inputFileNameLength=`expr $inputFile|wc -c`;
inputFileNameLength=`expr $inputFileNameLength - 1`;
echo "#| length of input file name   = "$inputFileNameLength >> $logFile;
fileName=${inputFile%.*};
echo "#| trimmed file name           = "$fileName >> $logFile;
fileNameLength=`expr $fileName|wc -c`;
fileNameLength=`expr $fileNameLength - 1`;
echo "#| length of trimmed file name = "$fileNameLength >> $logFile;
fileExtLength=`expr $inputFileNameLength - $fileNameLength - 1`;
echo "#| length of extension         = "$fileExtLength >> $logFile;
fileExt1=`echo $inputFile | sed -e "s/$fileName.//g"`;
echo "#| file extension              = "$fileExt1 >> $logFile;
fileExt2=$(echo "$fileExt1" | awk '{print tolower($0)}');	# convert to lowercase.
echo "#| file extension, lowercase   = "$fileExt1 >> $logFile;
# Check if inputfile is BAM format, then decompress if needed.
if [ "$fileExt2" = "bam" ]
then
	samtools view -h $projectDirectory$inputFile > $projectDirectory$fileName.sam
	newExt=".sam";
else
	newExt="."$fileExt1;
fi
echo "#| fileExt1 = "$fileExt1 >> $logFile;
echo "#| fileName = "$fileName >> $logFile;
echo "#| newExt   = "$newExt   >> $logFile;

tempOutput1=$projectDirectory"data_r1.b.fastq"
tempOutput2a=$projectDirectory"data_r2.b.fastq"
tempOutput2b=$projectDirectory"data_r2.c.fastq"

finalOutput1=$projectDirectory"data_r1.fastq";
finalOutput2=$projectDirectory"data_r2.fastq";

# Convert SAM to fastq using generalized commandline tools.
(
grep -v ^@ $projectDirectory$fileName$newExt | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > $tempOutput1;
) &

(
grep -v ^@ $projectDirectory$fileName$newExt | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > $tempOutput2a;
) &

wait;

# The resulting files do not match those from picard-tools in the following ways:
#--------------------------------------------------------------------------------
echo "#| Correct some formatting errors in the resulting FASTQ files." >> $logFile;
echo "#|     forward read file: append '/1' at end of each id line string." >> $logFile;
echo "#|     reverse read file: append '/2' at end of each id line string." >> $logFile;
echo "#|                        reverse-complement each sequence line string." >> $logFile;
echo "#|                        reverse each quality line string." >> $logFile;
## Correct the first-end read file.
(
# Add '/1' to each ID string.
sed 's/$/\/1/;n;n;n' $tempOutput1 > $finalOutput1;
) &

## Correct the second-end read file.
(
# Add '/2' to each ID string.
sed 's/$/\/2/;n;n;n' $tempOutput2a > $tempOutput2b;

# Reverse compliment fastq entries.
$seqtkDirectory"seqtk" seq -r $tempOutput2b > $finalOutput2;
) &

wait;


#===================================================================================================================================
# Delete the intermediate files which we don't want cluttering the drive.
#-----------------------------------------------------------------------------------------------------------------------------------
echo "#| Delete intermediate FASTQ files." >> $logFile;
rm $tempOutput1;
rm $tempOutput2a;
rm $tempOutput2b;

echo "#|---- sam2fastw.sh ---- end." >> $logFile;
