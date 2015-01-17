#!/bi
#n/bash -e
#
# genome.install_6.sh
#
set -e;
## All created files will have permission 760
umask 007;

## define script file locations.
user=$1;
genome=$2;
main_dir=$(pwd)"/../";

#user="darren";
#genome="Candida_albicans_SC5314_verA21-s02-m09-r07";
#main_dir="/heap/hapmap/bermanlab/";

reflocation=$main_dir"users/"$user"/genomes/"$genome"/";				# Directory where FASTA file is kept.
FASTA=`sed -n 1,1'p' $reflocation"reference.txt"`;						# Name of FASTA file.
FASTAname=$(echo $FASTA | sed 's/.fasta//g');							# name of genome file, without file type.
ddRADseq_FASTA=$FASTAname".MfeI_MboI.fasta";							# Name of digested reference for ddRADseq analysis.
RNAseq_FASTA=$FASTAname".expression.fasta";								# Name of digested reference for expression analysis.
standard_bin_FASTA=$FASTAname".standard_bins.fasta";					# Name of reference genome broken up into standard bins.

logName=$reflocation"process_log.txt";
condensedLog=$reflocation"condensed_log.txt";
#chmod 0755 $logName;
echo "\n\nRunning 'sh/genome.install_6.sh'" >> $logName;
echo "\tInput to shell script:" >> $logName;
echo "\t\t\$1 (user)         = $1" >> $logName;
echo "\t\t\$2 (genome)       = $2" >> $logName;
echo "" >> $logName;
echo "\tImportant location variables in script:" >> $logName;
echo "\t\t\$logName          = "$logName >> $logName;
echo "\t\t\$reflocation      = "$reflocation >> $logName;
echo "\t\t\$FASTA            = "$FASTA >> $logName;
echo "\t\t\$FASTAname        = "$FASTAname >> $logName;
echo "\t\t\$ddRADseq_FASTA   = "$ddRADseq_FASTA >> $logName;
echo "\t\t\$RNAseq_FASTA     = "$RNAseq_FASTA >> $logName;
echo "" >> $logName;
echo "Setting up for processing." >> $condensedLog;

# load local installed program location variables.
. $main_dir/sh/local_installed_programs.sh;

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
	bowtie2-build $reflocation$FASTA $reflocation"bowtie_index";
else
	echo "\tBowtie index for genome '$genome' found" >> $logName;
fi

#echo "\t============================================================================================== 2" >> $logName;
#
## Check if BLAST database has been made for selected genome. Generate database if not found.
#if [ -e $reflocation$FASTA".nin" ]
#then
#	echo "\tBLAST database for genome '$genome' found." >> $logName;
#else
#	echo "Generating BLAST database for genome." >> $condensedLog;
#	echo "\tBLAST database for genome '$genome' not found: Regenerating database." >> $logName;
#	formatdb -i $reflocation$FASTA -p F -o T;
#fi

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
	java -jar $picardDirectory"CreateSequenceDictionary.jar" R=$reflocation$FASTA O=$reflocation$FASTAname".dict";
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

# echo "\tRepetitiveness calculations have been removed from pipeline." >> $logName;
# echo "\tThey are very time-consuming and have been found to have little utility at present." >> $logName;

## Check if repetitiveness analysis has been done for genome.
repetgenome=$reflocation$FASTAname".repetitiveness.txt";

if [ -e $repetgenome ]
then
	echo "\tRepetitiveness file for genome '$genome' found." >> $logName;
else
	echo "Calculating repetitiveness of FASTA file for genome. (slow)" >> $condensedLog;
	echo "\tRepetitiveness file not found for genome '$genome': Regenerating using MatLab." >> $logName;

	## Perform repetitiveness analysis on reference file for genome.
	outputName=$reflocation"processing1.m";
	echo "\tWriting MATLAB function file to perform processing step." >> $logName;
	echo "\t\toutputName = "$outputName >> $logName;

	echo "function [] = processing1()" > $outputName;
	echo "\tdiary('"$reflocation"matlab.repetitiveness_analysis.log');" >> $outputName;
	echo "\tcd "$main_dir"Matlab/genome_install;" >> $outputName;
	echo "\trepetitiveness_1('$user','$genome');" >> $outputName;
	echo "\texit;" >> $outputName;
	echo "end" >> $outputName;

	echo "\t|\tfunction [] = processing1()" >> $logName;
	echo "\t|\t\tdiary('"$reflocation"matlab.repetitiveness_analysis.log');" >> $logName;
	echo "\t|\t\tcd "$main_dir"Matlab/genome_install;" >> $logName;
	echo "\t|\t\trepetitiveness_1('$user','$genome');" >> $logName;
	echo "\t|\t\texit;" >> $logName;
	echo "\t|\tend" >> $logName;

	echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
	matlab -nosplash -r "run "$outputName";";
	sed 's/^/\t\t\t|/;' $reflocation"matlab.repetitiveness_analysis.log" >> $logName;
fi

echo "\t============================================================================================== 6" >> $logName;

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
    echo "\texit;" >> $outputName;
    echo "end" >> $outputName;

    echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
    matlab -nosplash -r "run "$outputName";";
    sed 's/^/\t\t\t|/;' $reflocation"matlab.standard_fragmentation_of_reference.log" >> $logName;
fi

echo "\t----------------------------------------------------------------------------------------------" >> $logName;

if [ -e $reflocation$ddRADseq_FASTA ]
then
	echo "\tSimulated digest of genome already complete." >> $logName;
else
	echo "Performing simulated restriction digest of genome." >> $condensedLog;
	echo "\tSimulated restriction digest of genome being performed." >> $logName;

	## Perform simulated digest of genome.
	outputName=$reflocation"processing2.m";
	echo "\t\tWriting MATLAB function file to perform processing step." >> $logName;
	echo "\t\toutputName = "$outputName >> $logName;

	echo "function [] = processing2()" > $outputName;
	echo "\tdiary('"$reflocation"matlab.simulated_digest_of_reference.log');" >> $outputName;
	echo "\tcd "$main_dir"Matlab/genome_install;" >> $outputName;
	echo "\tgenome_process_for_RADseq_1('$user','$genome');" >> $outputName;
	echo "\texit;" >> $outputName;
	echo "end" >> $outputName;

	echo "\t|\tfunction [] = processing2()" >> $logName;
    echo "\t|\t\tdiary('"$reflocation"matlab.simulated_digest_of_reference.log');" >> $logName;
    echo "\t|\t\tcd "$main_dir"Matlab/genome_install;" >> $logName;
    echo "\t|\t\tgenome_process_for_RADseq_1('$user','$genome');" >> $logName;
    echo "\t|\t\texit;" >> $logName;
    echo "\t|\tend" >> $logName;

	echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
	matlab -nosplash -r "run "$outputName";";
	sed 's/^/\t\t\t|/;' $reflocation"matlab.simulated_digest_of_reference.log" >> $logName;
fi

echo "\t============================================================================================== 7" >> $logName;

inputFile=$reflocation"chromosome_features.txt";
outputFile=$reflocation"chromosome_features_2.txt";
if [ -e $inputFile ]
then
	echo "Simplifying and sorting chromosome_features file." >> $condensedLog;
	## Simplifying and sorting chromosome_features file.
	echo "\n\tSimplifying and sorting chromosome features file." >> $logName;
	echo "\n\t\tfeatures file = "$reflocation"chromosome_features.txt" >> $logName;
	echo "" > $outputFile;
	$python_exec $main_dir"py/genome/chromosome_features.simplify.py" $user $genome $main_dir $logName >> $outputFile;
else
	echo "\n\tChromosome features file not available." >> $logName;
fi

echo "\t----------------------------------------------------------------------------------------------" >> $logName;

if [ -e $reflocation$RNAseq_FASTA ]
then
	echo "\tExpression digest of genome already complete." >> $logName;
else
	if [ -e $reflocation"expression.txt" ]
	then
		echo "Performing simulated digest of genome into expression units." >> $condensedLog;
		echo "\tExpression digest of genome being performed." >> $logName;

		## Perform expression digest of genome.
		outputName=$reflocation"processing3.m";
		echo "\t\tWriting MATLAB function file to perform processing step." >> $logName;
		echo "\t\toutputName = "$outputName >> $logName;

		echo "function [] = processing3()" > $outputName;
		echo "\tdiary('"$reflocation"matlab.expression_digest_of_reference.log');" >> $outputName;
		echo "\tcd "$main_dir"Matlab/genome_install;" >> $outputName;
		echo "\tgenome_process_for_expression_1('$user','$genome');" >> $outputName;
		echo "\texit;" >> $outputName;
		echo "end" >> $outputName;

		echo "\t|\tfunction [] = processing3()" >> $logName;
		echo "\t|\t\tdiary('"$reflocation"matlab.expression_digest_of_reference.log');" >> $logName;
		echo "\t|\t\tcd "$main_dir"Matlab/genome_install;" >> $logName;
		echo "\t|\t\tgenome_process_for_expression_1('$user','$genome');" >> $logName;
		echo "\t|\t\texit;" >> $logName;
		echo "\t|\tend" >> $logName;

		echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
		matlab -nosplash -r "run "$outputName";";
		sed 's/^/\t\t\t|/;' $reflocation"matlab.expression_digest_of_reference.log" >> $logName;
	fi
fi

echo "\t============================================================================================== 7" >> $logName;

## Reformat standard-bin fragmented FASTA file to have single-line entries for each sequence fragment.
echo "Reformatting standard genome fragments FASTA file." >> $condensedLog;
echo "\tReformatting digested FASTA file => single-line per sequence fragment." >> $logName;
sh $main_dir"sh/FASTA_reformat_1.sh" $reflocation$standard_bin_FASTA;

outputFile=$reflocation$FASTAname".GC_ratios.standard_bins.txt";
if [ -e $outputFile ]
then
	echo "\n\tGC-ratios per standard-bin fragment has been calculated." >> $logName
else
	echo "Calculating GC ratios for genome standard-bin fragments." >> $condensedLog;
	## Calculating GC_ratio of standard bin fragments.
	echo "\n\tCalculating GC-ratios per each standard bin fragment." >> $logName;
	echo "\n\t\treflocation = "$reflocation >> $logName;
	echo "" > $outputFile;
	$python_exec $main_dir"py/genome/genome_process_for_standard_bins.GC_bias_1.py" $reflocation $logName >> $outputFile;
fi

if [ -e $repetgenome ]
then
	# This depends on whole genome repetitiveness analysis done previously.
	outputFile=$reflocation$FASTAname".repetitiveness.standard_bins.txt";
	if [ -e $outputFile ]
	then
		echo "\n\trepetitiveness per standard-bin fragment has been calculated." >> $logName
	else
		echo "Calculating repetitiveness of genome standard-bin fragments." >> $condensedLog;
		## Calculating repetitiveness of standard bin fragments.
		echo "\tCalculating repetitiveness per each digestion fragment." >> $logName;
		inputFile=$reflocation$FASTAname".repetitiveness.txt";
		echo "" > $outputFile;
		$python_exec $main_dir"py/genome/genome_process_for_standard_bins.repetitiveness_1.py" $inputFile $reflocation $logName >> $outputFile;
	fi
fi

echo "\t----------------------------------------------------------------------------------------------" >> $logName;

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
	$python_exec $main_dir"py/genome/genome_process_for_RADseq.GC_bias_1.py" $reflocation $logName >> $outputFile;
fi

if [ -e $repetgenome ]
then
	# This depends on whole genome repetitiveness analysis done previously.
	outputFile=$reflocation$FASTAname".repetitiveness.MfeI_MboI.txt";
	if [ -e $outputFile ]
	then
		echo "\n\tRepetitiveness per digestion fragment has been calculated." >> $logName
	else
		echo "Calculating repetitiveness for digested genome fragments." >> $condensedLog;
		## Calculating repetitiveness of ddRADseq (MfeI & MboI) fragments.
		echo "\n\n\tCalculating repetitiveness per each digestion fragment." >> $logName;
		inputFile=$reflocation$FASTAname".repetitiveness.txt";
		echo "" > $outputFile;
		$python_exec $main_dir"py/genome/genome_process_for_RADseq.repetitiveness_1.py" $inputFile $reflocation $logName >> $outputFile;
	fi
fi

echo "\t----------------------------------------------------------------------------------------------" >> $logName;

## Reformat digested FASTA file to have single-line entries for each sequence fragment.
echo "Reformatting expression genome fragments FASTA file." >> $condensedLog;
echo "\tReformatting expression FASTA file => single-line per sequence fragment." >> $logName;
sh $main_dir"sh/FASTA_reformat_1.sh" $reflocation$RNAseq_FASTA;

outputFile=$reflocation$FASTAname".GC_ratios.expression.txt";
if [ -e $outputFile ]
then
	echo "\n\tGC-ratios per expression fragment has been calculated." >> $logName
else
	echo "Calculating GC ratios for expression fragments." >> $condensedLog;
	## Calculating GC_ratio  of RNAseq (expression) fragments.
	echo "\n\tCalculating GC-ratios per each expression fragment." >> $logName;
	echo "\n\t\treflocation = "$reflocation >> $logName;
	echo "" > $outputFile;
	$python_exec $main_dir"py/genome/genome_process_for_RNAseq.GC_bias_1.py" $reflocation $logName >> $outputFile;
fi

if [ -e $repetgenome ]
then
	# This depends on whole genome repetitiveness analysis done previously.
	outputFile=$reflocation$FASTAname".repetitiveness.expression.txt";
	if [ -e $outputFile ]
	then
		echo "\n\tRepetitiveness per expression fragment has been calculated." >> $logName
	else
		echo "Calculating repetitiveness for expression fragments." >> $condensedLog;
		## Calculating repetitiveness of RNAseq (expression) fragments.
		echo "\n\n\tCalculating repetitiveness per each expression fragment." >> $logName;
		inputFile=$reflocation$FASTAname".repetitiveness.txt";
		echo "" > $outputFile;
		$python_exec $main_dir"py/genome/genome_process_for_RNAseq.repetitiveness_1.py" $inputFile $reflocation $logName >> $outputFile;
	fi
fi

echo "\n\t============================================================================================== 8" >> $logName;


##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------
sh $main_dir"sh/cleaning_genome.sh" $user $genome $main_dir;
