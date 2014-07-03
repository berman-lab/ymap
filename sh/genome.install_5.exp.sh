#!/bi
#n/bash -e
#
# Initialization of genome into pipeline
# 	$1 : user
# 	$2 : genome
#	$3 : directory
set -e;
## All created files will have permission 760
umask 007;

## define script file locations.
user=$1;
genome=$2;
main_dir=$3;

#user="default";
#genome="Candida_albicans_SC5314__A21-s02-m08-r09";
#main_dir="/heap/hapmap/bermanlab/";

reflocation=$main_dir"users/"$user"/genomes/"$genome"/";				# Directory where FASTA file is kept.
FASTA=`sed -n 1,1'p' $reflocation"reference.txt"`;						# Name of FASTA file.
FASTAname=$(echo $FASTA | sed 's/.fasta//g');							# name of genome file, without file type.
ddRADseq_FASTA=$FASTAname".MfeI_MboI.fasta";							# Name of digested reference for ddRADseq analysis.
standard_bin_FASTA=$FASTAname".standard_bins.fasta";					# Name of reference genome broken up into standard bins.

logName=$reflocation"process_log.txt";
condensedLog=$reflocation"condensed_log.txt";
#chmod 0755 $logName;
echo "\n\nRunning 'sh/genome.install_5.sh'" >> $logName;

echo "\tInput to shell script:" >> $logName;
echo "\t\t\$1 (user)       = $1" >> $logName;
echo "\t\t\$2 (genome)     = $2" >> $logName;
echo "\t\t\$3 (directory)  = $3" >> $logName;
echo "" >> $logName;
echo "\tImportant location variables in script:" >> $logName;
echo "\t\t\$logName        = "$logName >> $logName;
echo "\t\t\$reflocation    = "$reflocation >> $logName;
echo "\t\t\$FASTA          = "$FASTA >> $logName;
echo "\t\t\$FASTAname      = "$FASTAname >> $logName;
echo "\t\t\$ddRADseq_FASTA = "$ddRADseq_FASTA >> $logName;
echo "" >> $logName;
echo "Setting up for processing." >> $condensedLog;

##============================================#
# Initialization of various programs below.   #
#============================================##

echo "\t============================================================================================== 1" >> $logName;

## Check is genome and index files for Bowtie are available: Exit if genome files not found; Generate index files if needed.
if [ ! -e $reflocation"bowtie_index.4.bt2" ]
then
	echo "Generating Bowtie2 index for genome." >> $condensedLog;
	echo "\tBowtie index for genome '$genome' not found: Reindexing genome." >> $logName;
	## Bowtie 2 commands:
	/home/dabbey/software/bowtie2-2.1.0/bowtie2-build $reflocation$FASTA $reflocation"bowtie_index";
else
	echo "\tBowtie index for genome '$genome' found" >> $logName;
fi

echo "\t============================================================================================== 2" >> $logName;

## Check if BLAST database has been made for selected genome. Generate database if not found.
if [ -e $reflocation$FASTA".nin" ]
then
	echo "\tBLAST database for genome '$genome' found." >> $logName;
else
	echo "Generating BLAST database for genome." >> $condensedLog;
	echo "\tBLAST database for genome '$genome' not found: Regenerating database." >> $logName;
	formatdb -i $reflocation$FASTA -p F -o T;
fi

echo "\t============================================================================================== 3" >> $logName;

## Check if GATK dictionary and index files have been made for selected genome. Generate these files if not found.
if [ -e $reflocation$FASTAname".dict" ]
then
	echo "\tFASTA dictionary file for genome '$genome' found." >> $logName;
else
	echo "Generating FASTA dictionary file for genome, step 1." >> $condensedLog;
	echo "\tFASTA dictionary file not found for genome '$genome': Regenerating using Picard-tools." >> $logName;
	echo "\tR="$reflocation$FASTA >> $logName;
	echo "\tO="$reflocation$FASTAname".dict" >> $logName;
	java -jar /home/dabbey/software/picard-tools-1/picard-tools-1.105/CreateSequenceDictionary.jar R=$reflocation$FASTA O=$reflocation$FASTAname".dict";
fi

echo "\t============================================================================================== 4" >> $logName;

## Check if Samtools FASTA index file is found.
if [ -e $reflocation$FASTA".fai" ]
then
	echo "\tFASTA index file for genome '$genome' found." >> $logName;
else
	echo "Generatiing FASTA dictionary file for genome, step2." >> $condensedLog;
	echo "\tFASTA index file not found for genome '$genome': Regenerating using SamTools." >> $logName;
	samtools faidx $reflocation$FASTA;
fi

echo "\t============================================================================================== 5" >> $logName;

echo "\tRepetitiveness calculations have been removed from pipeline." >> $logName;
echo "\tThey are very time-consuming and have been found to have little utility at present." >> $logName;

### Check if repetitiveness analysis has been done for genome.
#repetgenome=$reflocation$FASTAname".repetitiveness.txt";
#
#if [ -e $repetgenome ]
#then
#	echo "\tRepetitiveness file for genome '$genome' found." >> $logName;
#else
#	echo "\tRepetitiveness file not found for genome '$genome': Regenerating using MatLab." >> $logName;
#
#	## Perform repetitiveness analysis on reference file for genome.
#	outputName=$reflocation"processing1.m";
#	echo "\tWriting MATLAB function file to perform processing step." >> $logName;
#	echo "\t\toutputName = "$outputName >> $logName;
#
#	echo "function [] = processing1()" > $outputName;
#	echo "\tdiary('"$reflocation"matlab.repetitiveness_analysis.log');" >> $outputName;
#	echo "\tcd "$main_dir"Matlab/genome_install;" >> $outputName;
#	echo "\trepetitiveness_1('$user','$genome');" >> $outputName;
#	echo "end" >> $outputName;
#
#	echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
#	matlab -nosplash -r "run "$outputName";";
#	sed 's/^/\t\t\t|/;' $reflocation"matlab.repetitiveness_analysis.log" >> $logName;
#fi

echo "\t============================================================================================== 6" >> $logName;

if [ -e $reflocation$ddRADseq_FASTA ]
then
	echo "\tSimulated digest of genome already complete." >> $logName;
else
	echo "Performing simulated restriction digest of genome." >> $condensedLog;
	echo "\tSimulated digest of genome being performed." >> $logName;

	## Perform simulated digest of genome.
	outputName=$reflocation"processing2.m";
	echo "\t\tWriting MATLAB function file to perform processing step." >> $logName;
	echo "\t\toutputName = "$outputName >> $logName;
	echo "function [] = processing2()" > $outputName;
	echo "\tdiary('"$reflocation"matlab.simulated_digest_of_reference.log');" >> $outputName;
	echo "\tcd "$main_dir"Matlab/genome_install;" >> $outputName;
	echo "\tgenome_process_for_RADseq_1('$user','$genome');" >> $outputName;
	echo "end" >> $outputName;

	echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
	matlab -nosplash -r "run "$outputName";";
	sed 's/^/\t\t\t|/;' $reflocation"matlab.simulated_digest_of_reference.log" >> $logName;
fi

echo "\t---------------------------------------------------------------------------------------------- 7" >> $logName;

## Reformat digested FASTA file to have single-line entries for each sequence fragment.
echo "Reformatting digested genome fragments FASTA file." >> $condensedLog;
echo "\tReformatting digested FASTA file => single-line per sequence fragment." >> $logName;
sh $main_dir"sh/FASTA_reformat_1.sh" $reflocation$ddRADseq_FASTA;

outputFile=$reflocation$FASTAname".GC_ratios.MfeI_MboI.txt";
if [ -e $outputFile ]
then
	echo "\n\tGC-ratios per digestion fragment has been calculated." >> $logName
else
	echo "Calculating GC ratios for digested genome fragments." >> $condensedLog;
	## Calculating GC_ratio  of ddRADseq (MfeI & MboI) fragments.
	echo "\n\tCalculating GC-ratios per each restriction digestion fragment." >> $logName;
	echo "\n\t\treflocation = "$reflocation >> $logName;
	echo "" > $outputFile;
	python $main_dir"py/genome/genome_process_for_RADseq.GC_bias_1.py" $reflocation $logName >> $outputFile;
fi

#outputFile=$reflocation$FASTAname".repetitiveness.MfeI_MboI.txt";
#if [ -e $outputFile ]
#then
#	echo "\n\tRepetitiveness per digestion fragment has been calculated." >> $logName
#else
#	## Calculating repetitiveness of ddRADseq (MfeI & MboI) fragments.   This depends on whole genome repetitiveness analysis.
#	echo "\n\n\tCalculating repetitiveness per each digestion fragment." >> $logName;
#	inputFile=$reflocation$FASTAname".repetitiveness.txt";
#	echo "" > $outputFile;
#	python $main_dir"py/genome/genome_process_for_RADseq.repetitiveness_1.py" $inputFile $reflocation $logName >> $outputFile;
#fi

echo "\n\t============================================================================================== 8" >> $logName;

if [ -e $reflocation$standard_bin_FASTA ]
then
	echo "\tGenome already fragmented into standard bins." >> $logName;
else
	echo "Performing standard-bin fragmentation of genome." >> $condensedLog;
	echo "\tGenome being fragmentated into standard bins." >> $logName;

    ## Perform sim.
    outputName=$reflocation"processing3.m";
	echo "\t\tWriting MATLAB function file to perform processing step." >> $logName;
    echo "\t\toutputName = "$outputName >> $logName;

    echo "function [] = processing3()" > $outputName;
    echo "\tdiary('"$reflocation"matlab.standard_fragmentation_of_reference.log');" >> $outputName;
    echo "\tcd "$main_dir"Matlab/genome_install;" >> $outputName;
    echo "\tgenome_process_for_standard_bins_1('$user','$genome');" >> $outputName;
    echo "end" >> $outputName;

	echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
    matlab -nosplash -r "run "$outputName";";
	sed 's/^/\t\t\t|/;' $reflocation"matlab.standard_fragmentation_of_reference.log" >> $logName;
fi

echo "\t---------------------------------------------------------------------------------------------- 9" >> $logName;

	## Reformat standard-bin fragmented FASTA file to have single-line entries for each sequence fragment.
	echo "\tReformatting digested FASTA file => single-line per sequence fragment." >> $logName;
	sh $main_dir"sh/FASTA_reformat_1.sh" $reflocation$standard_bin_FASTA;

outputFile=$reflocation$FASTAname".GC_ratios.standard_bins.txt";
if [ -e $outputFile ]
then
	echo "\n\tGC-ratios per standard bin fragment has been calculated." >> $logName
else
	echo "Calculating GC ratios for genome fragments." >> $condensedLog;
	## Calculating GC_ratio of standard bin fragments.
	echo "\n\tCalculating GC-ratios per each standard bin fragment." >> $logName;
	echo "\n\t\treflocation = "$reflocation >> $logName;
	echo "" > $outputFile;
	python $main_dir"py/genome/genome_process_for_standard_bins.GC_bias_1.py" $reflocation $logName >> $outputFile;
fi

outputFile=$reflocation$FASTAname".repetitiveness.standard_bins.txt";
#if [ -e $outputFile ]
#then
#	echo "\n\trepetitiveness per standard bin fragment has been calculated." >> $logName
#else
#	## Calculating repetitiveness of standard bin fragments.   This depends on whole genome repetitiveness analysis.
#	echo "\tCalculating repetitiveness per each digestion fragment." >> $logName;
#	inputFile=$reflocation$FASTAname".repetitiveness.txt";
#	echo "" > $outputFile;
#	python $main_dir"py/genome/genome_process_for_standard_bins.repetitiveness_1.py" $inputFile $reflocation $logName >> $outputFile;
#fi

echo "\n\t============================================================================================== 10" >> $logName;

	echo "Cleaning up intermediate files, step 1." >> $condensedLog;
	echo "\tDeleting unneeded intermediate files." >> $logName;
	rm $reflocation"processing2.m";
	rm $reflocation"processing3.m";
	rm $reflocation"matlab.simulated_digest_of_reference.log";
	rm $reflocation"matlab.standard_fragmentation_of_reference.log";

	nameString1=`cat $reflocation"name.txt"`;
	echo "\tName string for genome = '"$nameString1"'" >> $logName;

	echo "\tOpening up permissions, so genome can be transferred to the 'default' system user by any admin." >> $logName;
	echo "Cleaning up intermediate files, step 2." >> $condensedLog;
	chmod 0755 $reflocation"annotations.txt";
	chmod 0755 $reflocation"centromere_locations.txt";
	chmod 0755 $reflocation"chromosome_sizes.txt";
	chmod 0755 $reflocation"figure_definitions.txt";
	chmod 0755 $reflocation"name.txt";
	chmod 0755 $reflocation"working.txt";
	chmod 0755 $reflocation"reference.txt";
	chmod 0755 $reflocation"citation.txt";
	chmod 0755 $reflocation"ploidy.txt";

	echo "Cleaning up intermediate files, step 3." >> $condensedLog;
	chmod 0755 $reflocation"bowtie_index.1.bt2";
	chmod 0755 $reflocation"bowtie_index.2.bt2";
	chmod 0755 $reflocation"bowtie_index.3.bt2";
	chmod 0755 $reflocation"bowtie_index.4.bt2";
	chmod 0755 $reflocation"bowtie_index.rev.1.bt2";
	chmod 0755 $reflocation"bowtie_index.rev.2.bt2";

	echo "Cleaning up intermediate files, step 4." >> $condensedLog;
	chmod 0755 $reflocation$FASTAname".dict";
	chmod 0755 $reflocation$FASTA;
	chmod 0755 $reflocation$FASTA".fai";
	chmod 0755 $reflocation$FASTA".nhr";
	chmod 0755 $reflocation$FASTA".nin";
	chmod 0755 $reflocation$FASTA".nsd";
	chmod 0755 $reflocation$FASTA".nsi";
	chmod 0755 $reflocation$FASTA".nsq";
	chmod 0755 $reflocation$ddRADseq_FASTA;
	chmod 0755 $reflocation$standard_bin_FASTA;

	echo "Cleaning up intermediate files, step 5." >> $condensedLog;
	chmod 0755 $reflocation$FASTAname".seq.mat";
	chmod 0755 $reflocation$FASTAname".seqRevCom.mat";

	echo "Concluding analysis." >> $condensedLog;
	## Output a simple text file to tell the pipeline system that the genome installation has completed.
	echo "\tGenerating 'complete.txt' file to let pipeline know installation of genome has completed." >> $logName;
	echo "complete" > $main_dir"users/"$user"/genomes/"$genome"/complete.txt";
	chmod 0755 $reflocation"complete.txt";

	chmod 0755 $reflocation"condensed_log.txt";

echo "\t============================================================================================== 11" >> $logName;
