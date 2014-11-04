# Accept input.
user=$1;
project=$2;
main_dir=$3;
inputFile=$4;

# Define directories.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
picardDirectory="/home/dabbey/software/picard-tools-1/picard-tools-1.105/";

# Run picard-tools SamToFastq.
java -Xmx2g -jar $picardDirectory"SamToFastq.jar" INPUT=$projectDirectory$inputFile FASTQ=$projectDirectory"data_r1.fastq" SECOND_END_FASTQ=$projectDirectory"data_r2.fastq";
