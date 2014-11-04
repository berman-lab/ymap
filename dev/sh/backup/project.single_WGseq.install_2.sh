#!/bin/bash -e
#
# Initialization of genome into pipeline
#   $1 : user
#   $2 : project
#   $3 : base directory
set -e;
## All created files will have permission 760
umask 007;

## define script file locations.
user=$1;
project=$2;
main_dir=$3;


##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------

# Define project directory.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";

# Setup process_log.txt file.
logName=$projectDirectory"process_log.txt";
chmod 0755 $logName;
echo "Log file initiated during processing of 'sh/project.single_WGseq.install_2.sh'\n" >> $logName;

echo "#============================================================================== 1\n" >> $logName;

echo "projectDirectory = '"$projectDirectory"'\n"" >> $logName;

# Get genome name and hapmap use from project's "genome.txt" file.
genome=$(head -n 1 $projectDirectory"genome.txt");
hapmapBoolean=$(tail -n 1 $projectDirectory"genome.txt");
echo "Genome = '"$genome"'\n" >> $logName;
if [ $genome -eq $hapmapBoolean ]
then
	# only one line, so hapmap is not available for the genome being examined.
	hapmapBoolean="F";
	echo "No hapmap is available.\n" >> $logName;
else
	# file has two different lines, so genome being examined has hapmap data.
	if [ $hapmapBoolean -eq "Y" ]
	then
		echo "A hapmap is available and is in use.\n" >> $logName;
	else
		echo "A hapmap is available, but is not in use.\n" >> $logName;
	fi
fi

# Determine location of genome being used.
if [ -d $main_dir"users/"$user"/genomes/"$genome"/" ]
then
	genomeDirectory=$main_dir"users/"$user"/genomes/"$genome"/";
elif [ -d $main_dir"users/default/genomes/"$genome"/" ]
	genomeDirectory=$main_dir"users/default/genomes/"$genome"/";
fi
echo "genomeDirectory = '"$genomeDirectory"'\n";" >> $logName;

# Get reference FASTA file name from "reference.txt";
genomeFASTA=$(head -n 1 $genomeDirectory"reference.txt");

# Get data file name from "datafiles.txt";
datafile=$(head -n 1 $projectDirectory"datafiles.txt");

# Define temporary directory for FASTQC files.
fastqcDirectory=$projectDirectory"fastqc_temp/";

echo "#============================================================================== 2\n" >> $logName;

##==============================================================================
## Initial processing of single-WGseq dataset.
##------------------------------------------------------------------------------

# Align fastq against genome.
echo "[[=- Align with Bowtie -=]]\n" >> $logName;
threads=8;

echo "\tBowtie : single-end reads aligning into SAM file.\n" >> $logName;
## Bowtie 2 command for single reads:
bowtie2 --very-sensitive -p $threads $genomeDirectory"bowtie_index_"$genome -U $datafile -S $projectDirectory"data.sam";
	# -S : SAM output mode.
	# -p : number of threads to use.
	# -1 : dataset.
    # --very-sensitive : a default set of configurations.
echo "\tBowtie : single-end reads aligned into SAM file.\n" >> $logName;

echo "\tSamtools : converting Bowtie-SAM into compressed format (BAM) file." >> $logName;
samtools view -bT $genomeDirectory$genomeFASTA $projectDirectory"data.sam" > $projectDirectory"data.temp.bam";
rm $projectDirectory"data.sam";
echo "\tSamtools : Bowtie-SAM converted into compressed format (BAM) file.\n" >> $logName;

echo "\tPicard : Adding headers to Bowtie-BAM file.\n" >> $logName;
java -Xmx2g -jar /soft/picard-tools/1.68/AddOrReplaceReadGroups.jar INPUT=$projectDirectory"data.temp.bam" OUTPUT=$projectDirectory"data.bam" RGID=1 RGLB=1 RGPL=ILLUMINA RGPU=1 RGSM=SM;
rm $projectDirectory"data.temp.bam";
echo "\tPicard : Headers added to Bowtie-BAM file.\n" >> $logName;

echo "[[=- Sorting/Indexing BAM files -=]]\n" >> $logName;

echo "\tSamtools : Bowtie-BAM sorting, deduping, & indexing.\n" >> $logName;
samtools sort $projectDirectory"data.bam" $projectDirectory"data_sorted";
samtools index $projectDirectory"data_sorted.bam";
echo "\tSamtools : Bowtie-BAM sorted & indexed.\n" >> $logName;

echo "[[=- GATK analysis, indel-realignment -=]]\n" >> $logName;

inputFile=$projectDirectory"data_sorted.bam";
outputFile1=$projectDirectory"data_forIndelRealigner.intervals";
outputFile2=$projectDirectory"data_indelRealigned.bam";
#outputFile3=$projectDirectory"data_SNPs_GATK.vcf";
GATKreference=$genomeDirectory$genomeFASTA;

# 'LENIENT' should allow GATK to ignore problem reads.
GATKoptions="-S LENIENT -filterMBQ";

## FASTQC : determine read quality format in use.
echo "\tFASTQC : Read quality coding is required.\n" >> $logName;
##-------------------------------------------------------------------------------------------------------
## FASTQC analysis of FASTQ input files to determine read quality coding.
##      This is needed for determining GATK command option to deal with alternate coding.
##-------------------------------------------------------------------------------------------------------
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
echo "\t\tWorking directory : '"$fastqcDirectory"'\n" >> $logName;

# -o : points to standardized ouptut temp directory.
fastqc -o $fastqcDirectory $projectDirectory$datafile;

# get file extension.
fileExt=${datafile#*.};

# generate file name with extension removed and '_fastqc' appended, to match FASTQC output files.
FASTQC_temp_file="${datafile/.$fileExt/_fastqc}";

# parse read quality encoding format from FASTQC output.
qualityCodingLine=`sed -n 6,6'p' $fastqcDirectory$FASTQC_temp_file"/fastqc_data.txt"`;
qualityCoding=${qualityCodingLine:9};

# cleanup extraneous FASTQC files.
rm $fastqcDirectory$FASTQC_temp_file".zip";
rm -rf $$fastqcDirectory$FASTQC_temp_file;

# Trim starting and ending whitespace from quality coding string.
trimmedQualityCoding="${qualityCoding#"${qualityCoding%%[![:space:]]*}"}";                # remove leading whitespace characters
trimmedQualityCoding="${trimmedQualityCoding%"${trimmedQualityCoding##*[![:space:]]}"}";  # remove trailing whitespace characters
echo "\t\tQuality coding method = '"$trimmedQualityCoding"'.\n" >> $logName;

# Expected output:
#       "Illumina 1.5"          : starting with 64.
#       "Sanger / Illumina 1.9" : normal, starting with 33.   [MiSeq]

# If quality coding string matches a type known to start with 64, update GATK option string.
if [ "$trimmedQualityCoding" == "Illumina 1.5" ]
then
	GATKoptions=$GATKoptions" -fixMisencodedQuals";
fi
echo "\t\tGATK options          = '"$GATKoptions"'.\n" >> $logName;

echo "\tGATK : preparing for IndelRealignment.\n" >> $logName;
gatk -nt $threads -T RealignerTargetCreator -I $inputFile -R $GATKreference -o $outputFile1 $GATKoptions;
echo "\tGATK : prepared for IndelRealignment.\n" >> $logName;

echo "\tGATK : performing IndelRealignment." >> $logName;
gatk -T IndelRealigner -I $inputFile -R $GATKreference -targetIntervals $outputFile1 -o $outputFile2 $GATKoptions;
echo "\tGATK : performed IndelRealignment." >> $logName;

echo "#============================================================================== 3\n" >> $logName;

echo "[[=- In-house SNP/CNV/INDEL analysis -=]]\n" >> $logName;
usedFile=$projectDirectory"data_indelRealigned.bam";

echo "\tSamtools : Generating pileup.   (for SNP/CNV/INDEL analysis)\n" >> $logName;
samtools mpileup -f $genomeDirectory$genomeFASTA $usedFile | awk '{print $1 " " $2 " " $3 " " $4 " " $5}' > $projectDirectory"data.pileup";

# Run the in-house section before GATK and use the following line if GATK crashes.
#samtools mpileup -f $genomeDirectory$genomeFASTA $projectDirectory"data_cleaned.bam" | awk '{print $1 " " $2 " " $3 " " $4 " " $5}' > $projectDirectory"data.pileup";
# If GATK is working, unhide the next line.
# samtools mpileup -f $genomeDirectory$genomeFASTA $projectDirectory"data_indelRealigned.bam" | awk '{print $1 " " $2 " " $3 " " $4 " " $5}' > $projectDirectory"data.pileup";
echo "\tSamtools : Pileup generated.\n" >> $logName;

( echo "\tPython : Processing pileup for SNPs.\n" >> $logName;
python $main_dir"py/counts_SNPs_v5.py" $projectDirectory"data.pileup" > $projectDirectory"putative_SNPs_v4.txt";
echo "\tPython : Pileup processed for SNPs.\n" >> $logName; ) &

( echo "\tPython : Processing pileup for CNVs.\n" >> $logName;
python $main_dir"py/counts_CNVs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_CNVs_v1.txt";
echo "\tPython : Pileup processed for CNVs.\n" >> $logName; ) &

( echo "\tPython : Processing pileup for INDELs.\n" >> $logName;
python $main_dir"py/counts_INDELs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_INDELS_v1.txt";
echo "\tPython : Pileup processed for INDELs.\n" >> $logName; ) &

( echo "\tPython : Processing pileup for SNP-CNV.\n" >> $logName;
python $main_dir"py/counts_CNVs-SNPs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"SNP_CNV_v1.txt";
echo "\tPython : Pileup processed for SNP-CNV.\n" >> $logName; ) &

wait;

## Generate "complete.txt" to indicate processing has completed normally.
completeFile=$projectDirectory"complete.txt";
echo "complete" > $completeFile;
