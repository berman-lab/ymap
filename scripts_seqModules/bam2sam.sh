#===================================================================================================================================
# Convert Bam file to Sam file for introduction into the sequence analysis pipeline.
#-----------------------------------------------------------------------------------------------------------------------------------
user=$1;
project=$2;
inputFile=$3;
main_dir=$(pwd)"/";

projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logFile=$projectDirectory"process_log.txt";
echo "#|---- bam2sam.sh ---- begin." >> $logFile;

outputFile="data.sam"
echo "#| input file  : "$inputFile >> $logFile;
echo "#| output file : "$outputFile >> $logFile;

#===================================================================================================================================
# Use SAMtools to convert Bam to Sam.
#-----------------------------------------------------------------------------------------------------------------------------------

cd $projectDirectory;

echo "#| indexing BAM file." >> $logFile;
samtools index $inputFile;
echo "#| converting to SAM file." >> $logFile;
samtools view -h $inputFile > $outputFile;
echo "#|---- bam2sam.sh ---- end." >> $logFile;
