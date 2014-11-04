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

seqtkDirectory="/home/dabbey/software/Seqtk/";

# #===================================================================================================================================
# # Run picard-tools SamToFastq : This method requires large amounts of memory and will crash when it runs out of memory.
# #-----------------------------------------------------------------------------------------------------------------------------------
# picardDirectory="/home/dabbey/software/picard-tools-1/picard-tools-1.105/";
# java -Xmx16g -jar $picardDirectory"SamToFastq.jar" INPUT=$projectDirectory$inputFile FASTQ=$projectDirectory"data_r1.a.fastq" SECOND_END_FASTQ=$projectDirectory"data_r2.a.fastq";


#===================================================================================================================================
# Alternate method using generalized commandline tools with low memory footprint.
#-----------------------------------------------------------------------------------------------------------------------------------
echo "#| Translating SAM/BAM file to FASTQ format." >> $logFile;
fileExt1=${inputFile#*.};									# get file extension.
fileName=${inputFile%.*};									# trim extension off of file.
fileExt2=$(echo "$fileExt1" | awk '{print tolower($0)}');	# convert to lowercase.
# Check if inputfile is BAM format, then decompress if needed.
if [ "$fileExt2" = "bam" ]
then
	samtools view -h $projectDirectory$inputFile > $projectDirectory$fileName.sam
	newExt=".sam";
else
	newExt="."$fileExt1;
fi

Boutfile1=$projectDirectory"data_r1.b.fastq"
Boutfile2=$projectDirectory"data_r2.b.fastq"
Coutfile2=$projectDirectory"data_r2.c.fastq"

# Convert SAM to fastq using generalized commandline tools.
(
grep -v ^@ $projectDirectory$fileName$newExt | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > $Boutfile1;
) &

(
grep -v ^@ $projectDirectory$fileName$newExt | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > $Boutfile2;
) &

wait;

# The resulting files do not match those form picard-tools in the following ways:
#--------------------------------------------------------------------------------
# data_r1.b.fastq:
#     Append "/1" to each id string.
# data_r2.b.fastq:
#     Append "/2" to each id string.
#     Reverse-complement each sequence string.
#     Reverse each quality string.

output1=$projectDirectory"data_r1.fastq";
output2=$projectDirectory"data_r2.fastq";

echo "#| Correct some formatting errors in the resulting FASTQ files." >> $logFile;
echo "#|     forward read file: append '/1' at end of each id line string." >> $logFile;
echo "#|     reverse read file: append '/2' at end of each id line string." >> $logFile;
echo "#|                        reverse-complement each sequence line string." >> $logFile;
echo "#|                        reverse each quality line string." >> $logFile;
## Correct the first-end read file.
(
exec 4< $Boutfile1;			# open extra file descriptor for input.
while read S_id_line <&4;	# read lines from input file in blocks of four, then adjust them to match the expected format.
do
	read seq_line <&4;
	read Q_id_line <&4;
	read qual_line <&4;
	echo $S_id_line"/1" >> $output1;
	echo $seq_line      >> $output1;
	echo $Q_id_line     >> $output1;
	echo $qual_line     >> $output1;
done
exec 4<&-;	# close the extra file descriptors.
) &

## Correct the second-end read file.
(
exec 5< $Boutfile2;
while read S_id_line <&5;
do
	read seq_line <&5;
	read Q_id_line <&5;
	read qual_line <&5;
	echo $S_id_line"/2" >> $Coutfile2;
	echo $seq_line      >> $Coutfile2;
	echo $Q_id_line     >> $Coutfile2;
	echo $qual_line     >> $Coutfile2;
done
exec 5<&-;

$seqtkDirectory"seqtk" seq -r $Coutfile2 > $output2;

) &

wait;


#===================================================================================================================================
# Delete the intermediate files which we don't want cluttering the drive.
#-----------------------------------------------------------------------------------------------------------------------------------
echo "#| Delete intermediate FASTQ files." >> $logFile;
rm $Boutfile1;
rm $Boutfile2;
rm $Coutfile2;

echo "#|---- sam2fastw.sh ---- end." >> $logFile;
