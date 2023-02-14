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

reflocation=$main_dir"users/"$user"/genomes/"$genome"/";				# Directory where FASTA file is kept.
FASTA=`sed -n 1,1'p' $reflocation"reference.txt"`;					# Name of FASTA file.
FASTAname=$(echo $FASTA | sed 's/\.fasta//g');						# Name of genome file, without file type.
FASTA2=$(echo $FASTA | sed 's/\.fasta/\.2\.fasta/g');					# Name of reformatted genome file, to single-line entries.
repetgenome=$reflocation$FASTAname".repetitiveness.txt";				# Name of repetitiveness profile for genome.
repetgenome_smoothed=$reflocation$FASTAname".repetitiveness_smoothed.txt";		# Name of Gaussian smoothed repetitiveness profile for genome.
standard_bin_FASTA=$reflocation$FASTAname".standard_bins.fasta";			# Name of reference genome broken up into standard bins.
ddRADseq_FASTA=$reflocation$FASTAname".MfeI_MboI.fasta";				# Name of digested reference for ddRADseq analysis.
logName=$reflocation"process_log.txt";
condensedLog=$reflocation"condensed_log.txt";

#chmod 0755 $logName;
echo "\n\nRunning 'scripts_genomes/genome.install_6.sh'" >> $logName;
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
echo "" >> $logName;
echo "Setting up for processing." >> $condensedLog;

# load local installed program location variables.
. $main_dir"local_installed_programs.sh";

##============================================#
# Initialization of various programs below.   #
#============================================##

echo "\n\t============================================================================================== 1" >> $logName;

## Check is genome and index files for Bowtie are available: Exit if genome files not found; Generate index files if needed.
if [ ! -e $reflocation"bowtie_index.4.bt2" ]
then
	echo "Generating Bowtie2 index for genome." >> $condensedLog;
	echo "\tBowtie index for genome '$genome' not found: Reindexing genome." >> $logName;
	## Bowtie 2 commands:
	$bowtie2Directory"bowtie2-build" $reflocation$FASTA $reflocation"bowtie_index";
else
	echo "\tBowtie index for genome '$genome' found" >> $logName;
fi

#echo "\n\t============================================================================================== 2" >> $logName;
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

echo "\n\t============================================================================================== 3" >> $logName;

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

echo "\n\t============================================================================================== 4" >> $logName;

## Check if Samtools FASTA index file is found.
if [ -e $reflocation$FASTA".fai" ]
then
	echo "\tFASTA index file for genome '$genome' found." >> $logName;
else
	echo "Generatiing FASTA dictionary file for genome, step2." >> $condensedLog;
	echo "\tFASTA index file not found for genome '$genome': Regenerating using SamTools." >> $logName;
	$samtools_exec faidx $reflocation$FASTA;
fi

echo "\n\t============================================================================================== 5" >> $logName;

## Generate version of FASTA genome file to have single-line entries.
echo "\tReformatting genome FASTA file into single-line entries." >> $logName;
cp $reflocation$FASTA $reflocation$FASTA2;
sh $main_dir"scripts_seqModules/FASTA_reformat_1.sh" $reflocation$FASTA2;

echo "\n\t============================================================================================== 5" >> $logName;

## Check if repetitiveness analysis has been done for genome.
if [ -e $repetgenome ]
then
	echo "\tRepetitiveness file for genome '$genome' found." >> $logName;
else
	echo "Calculating repetitiveness of FASTA file for genome." >> $condensedLog;
	echo "\tRepetitiveness file not found for genome '$genome': Regenerating using Python script." >> $logName;

	## Perform repetitiveness analysis on reference file for genome, then smooth the profile.
	echo "" > $repetgenome;
        $python_exec $main_dir"scripts_genomes/repetitiveness_1.py"      $user $genome $main_dir $logName     >> $repetgenome;
	echo "" > $repetgenome_smoothed;
	$python_numpy_exec $main_dir"scripts_genomes/repetitiveness_smooth.py" $user $genome $main_dir $logName 128 >> $repetgenome_smoothed;
	mv $repetgenome_smoothed $repetgenome;
fi

echo "\n\t============================================================================================== 6" >> $logName;

if [ -e $standard_bin_FASTA ]
then
	echo "\tGenome already fragmented into standard bins." >> $logName;
else
	echo "Performing standard-bin fragmentation of genome." >> $condensedLog;
	echo "\tGenome being fragmentated into standard bins." >> $logName;

	## Perform reference genome fragmentation.
	echo "" > $standard_bin_FASTA;
	$python_exec $main_dir"scripts_genomes/genome_process_for_standard_bins_1.py" $user $genome $main_dir $logName >> $standard_bin_FASTA;
fi

echo "\n\t----------------------------------------------------------------------------------------------" >> $logName;

if [ -e $ddRADseq_FASTA ]
then
	echo "\tSimulated restriction digest (MfeI & MboI) of genome already complete." >> $logName;
else
	echo "Performing simulated restriction digest (MfeI & MboI) of genome." >> $condensedLog;
	echo "\tSimulated restriction digest of genome being performed." >> $logName;

	## Perform simulated digest of genome.
	echo "" > $ddRADseq_FASTA;
	$python_exec $main_dir"scripts_genomes/genome_process_for_RADseq_1.py" $user $genome $main_dir $logName >> $ddRADseq_FASTA;
fi

echo "\n\t============================================================================================== 7" >> $logName;

inputFile=$reflocation"chromosome_features.txt";
outputFile=$reflocation"chromosome_features_2.txt";
if [ -e $inputFile ]
then
	echo "Simplifying and sorting chromosome_features file." >> $condensedLog;
	## Simplifying and sorting chromosome_features file.
	echo "\n\tSimplifying and sorting chromosome features file." >> $logName;
	echo "\n\t\tfeatures file = "$reflocation"chromosome_features.txt" >> $logName;
	echo "" > $outputFile;
	$python_exec $main_dir"scripts_genomes/chromosome_features.simplify.py" $user $genome $main_dir $logName >> $outputFile;
else
	echo "\n\tChromosome features file not available." >> $logName;
fi

echo "\n\t============================================================================================== 7" >> $logName;

## Reformat standard-bin fragmented FASTA file to have single-line entries for each sequence fragment.
echo "Reformatting standard genome fragments FASTA file." >> $condensedLog;
echo "\tReformatting digested FASTA file => single-line per sequence fragment." >> $logName;
sh $main_dir"scripts_seqModules/FASTA_reformat_1.sh" $standard_bin_FASTA;

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
	$python_exec $main_dir"scripts_genomes/genome_process_for_standard_bins.GC_bias_1.py" $user $genome $main_dir $logName >> $outputFile;
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
		echo "\n\tCalculating repetitiveness per each digestion fragment." >> $logName;
		inputFile=$reflocation$FASTAname".repetitiveness.txt";
		echo "" > $outputFile;
		$python_exec $main_dir"scripts_genomes/genome_process_for_standard_bins.repetitiveness_2.py" $user $genome $main_dir $logName >> $outputFile;
	fi
fi

echo "\n\t----------------------------------------------------------------------------------------------" >> $logName;

## Reformat digested FASTA file to have single-line entries for each sequence fragment.
echo "Reformatting digested genome fragments FASTA file." >> $condensedLog;
echo "\tReformatting digested FASTA file => single-line per sequence fragment." >> $logName;
sh $main_dir"scripts_seqModules/FASTA_reformat_1.sh" $ddRADseq_FASTA;

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
	$python_exec $main_dir"scripts_genomes/genome_process_for_RADseq.GC_bias_1.py" $user $genome $main_dir $logName >> $outputFile;
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
		echo "" > $outputFile;
		$python_exec $main_dir"scripts_genomes/genome_process_for_RADseq.repetitiveness_2.py" $user $genome $main_dir $logName >> $outputFile;
	fi
fi

echo "\n\t============================================================================================== 8" >> $logName;

##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------
sh $main_dir"scripts_genomes/cleaning_genome.sh" $user $genome $main_dir;
