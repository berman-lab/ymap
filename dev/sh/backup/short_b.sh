#!/bin/bash -e
#
# Initialization of genome into pipeline
#   $1 : user
#   $2 : project
#   $3 : base directory
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
#user=$1;
#project=$2;
#main_dir=$3;
user="darren";
project="test";
main_dir="/heap/hapmap/bermanlab/";

##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------

# locations of auxillary software for pipeline analysis.
bowtie2Directory="/home/dabbey/software/bowtie2-2.1.0/";
picardDirectory="/home/dabbey/software/picard-tools-1/picard-tools-1.105/";
fastqcDirectory="/home/dabbey/software/FastQC/";
gatkDirectory="/home/dabbey/software/gatk/GenomeAnalysisTK-2.8-1-g932cd3a/";
java7Directory="/home/dabbey/software/Java7/jdk1.7.0_51/jre/bin/";

# Define project directory.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";

# Setup process_log.txt file.
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";
#chmod 0755 $logName;
echo "#.............................................................................." >> $logName;
echo "Running 'sh/project.single_WGseq.install_3.sh'" >> $logName;
echo "Variables passed via command-line from 'php/project.single_WGseq.install_2.php' :" >> $logName;
echo "\tuser     = '"$user"'" >> $logName;
echo "\tproject  = '"$project"'" >> $logName;
echo "\tmain_dir = '"$main_dir"'" >> $logName;
echo "#============================================================================== 3" >> $logName;

echo "#=====================================#" >> $logName;
echo "# Setting up locations and variables. #" >> $logName;
echo "#=====================================#" >> $logName;

echo "\tprojectDirectory = '$projectDirectory'" >> $logName;
echo "Setting up for processing." >> $condensedLog;

# Get genome name and hapmap use from project's "genome.txt" file.
genome=$(head -n 1 $projectDirectory"genome.txt");
hapmapBoolean=$(tail -n 1 $projectDirectory"genome.txt");
lineCount=$(wc -l < $projectDirectory"genome.txt");
echo "\t'genome.txt' file entries." >> $logName;
echo "\t\tlineCount = '"$lineCount"'" >> $logName;
echo "\t\tgenome = '"$genome"'" >> $logName;
if [ $lineCount -eq 1 ]
then
	# only one line, so hapmap is not available for the genome being examined.
	hapmapBoolean="N";
	echo "\tNo hapmap is available." >> $logName;
else
	echo "\t\thapmapBoolean = '"$hapmapBoolean"'" >> $logName;
	# file has two different lines, so genome being examined has hapmap data.
	if [ "$hapmapBoolean" == "Y" ]
	then
		echo "\tA hapmap is available and is in use." >> $logName;
	else
		echo "\tA hapmap is available, but is not in use." >> $logName;
	fi
fi

# Determine location of genome being used.
if [ -d $main_dir"users/"$user"/genomes/"$genome"/" ]
then
	genomeDirectory=$main_dir"users/"$user"/genomes/"$genome"/";
	genomeUser=$user;
elif [ -d $main_dir"users/default/genomes/"$genome"/" ]
then
	genomeDirectory=$main_dir"users/default/genomes/"$genome"/";
	genomeUser="default";
fi
echo "\tgenomeDirectory = '"$genomeDirectory"'" >> $logName;

# Get reference FASTA file name from "reference.txt";
genomeFASTA=$(head -n 1 $genomeDirectory"reference.txt");
echo "\tgenomeFASTA = '"$genomeFASTA"'" >> $logName;

# Get data file name from "datafiles.txt";
datafile=$(head -n 1 $projectDirectory"datafiles.txt");
echo "\tdatafile = '"$datafile"'" >> $logName;

# Define temporary directory for FASTQC files.
fastqcTempDirectory=$projectDirectory"fastqc_temp/";
echo "\tfastqcTempDirectory = '"$fastqcTempDirectory"'" >> $logName;

# Get ploidy estimate from "ploidy.txt" in project directory.
ploidyEstimate=$(head -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyEstimate = '"$ploidyEstimate"'" >> $logName;

# Get ploidy baseline from "ploidy.txt" in project directory.
ploidyBase=$(tail -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyBase = '"$ploidyBase"'" >> $logName;

# Get parent name from "parent.txt" in project directory.
projectParent=$(head -n 1 $projectDirectory"parent.txt");
echo "\tparentProject = '"$projectParent"'" >> $logName;

echo "#============================================================================== 2" >> $logName;

##==============================================================================
## Initial processing of single-WGseq dataset.
##------------------------------------------------------------------------------

##==============================================================================
## Perform CGH analysis, with GC-correction, on dataset.
##------------------------------------------------------------------------------

echo "#==========================#" >> $logName;
echo "# CGH analysis of dataset. #" >> $logName;
echo "#==========================#" >> $logName;
echo "Analyzing and mapping CNVs." >> $condensedLog;

echo "\t\tGenerating MATLAB script to perform CNV analysis of dataset, with GC-correction." >> $logName;
outputName=$projectDirectory"processing1.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing1()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.CNV_and_GCbias.log');" >> $outputName;
echo "\tcd "$main_dir"Matlab/WGseq;" >> $outputName;
echo "\tanalyze_CNVs_1('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t\tCalling MATLAB." >> $logName;
echo "================================================================================================";
echo "== CGH analysis ================================================================================";
echo "================================================================================================";
matlab -nosplash -r "run "$outputName"; exit;";
#chmod 0755 $outputName;
#chmod 0755 $projectDirectory"matlab.CNV_and_GCbias.log";
echo "\t\tMATLAB log from CNV analysis." >> $logName;
sed 's/^/\t\t\t|/;' $projectDirectory"matlab.CNV_and_GCbias.log" >> $logName;
rm $outputName;
rm $projectDirectory"matlab.CNV_and_GCbias.log";

##==============================================================================
## Perform ChARM analysis of dataset.
##------------------------------------------------------------------------------

echo "#============================#" >> $logName;
echo "# ChARM analysis of dataset. #" >> $logName;
echo "#============================#" >> $logName;
echo "Analyzing CNV breakpoints." >> $condensedLog;

echo "\t\tGenerating MATLAB script to perform ChARM analysis of dataset." >> $logName;
outputName=$projectDirectory"processing2.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing2()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.ChARM.log');" >> $outputName;
echo "\tcd "$main_dir"Matlab/ChARM;" >> $outputName;
echo "\tChARM_v4('$project','$user','$genome','$genomeUser','$main_dir');" >> $outputName;
echo "end" >> $outputName;

echo "\t\tCalling MATLAB." >> $logName;
echo "================================================================================================";
echo "== ChARM analysis ==============================================================================";
echo "================================================================================================";
matlab -nosplash -r "run "$outputName"; exit;";
#chmod 0755 $outputName;
#chmod 0755 $projectDirectory"matlab.ChARM.log";
echo "\t\tMATLAB log from ChARM analysis." >> $logName;
sed 's/^/\t\t\t|/;' $projectDirectory"matlab.ChARM.log" >> $logName;
rm $outputName;
rm $projectDirectory"matlab.ChARM.log";

##==============================================================================
## Perform SNP analysis on dataset.
##------------------------------------------------------------------------------

echo "#==========================#" >> $logName;
echo "# SNP analysis of dataset. #" >> $logName;
echo "#==========================#" >> $logName;

echo "Analyzing and mapping SNPs. " >> $condensedLog;

echo "\t\tGenerating MATLAB script to perform SNP analysis of dataset." >> $logName;
outputName=$projectDirectory"processing3.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing3()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.SNP_analysis.log');" >> $outputName;
echo "\tcd "$main_dir"Matlab/WGseq;" >> $outputName;
echo "\tanalyze_SNPs('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t\tCalling MATLAB." >> $logName;
echo "================================================================================================";
echo "== SNP analysis ================================================================================";
echo "================================================================================================";
matlab -nosplash -r "run "$outputName"; exit;";
#chmod 0755 $outputName;
#chmod 0755 $projectDirectory"matlab.SNP_analysis.log";
echo "\t\tMATLAB log from SNP analysis." >> $logName;
sed 's/^/\t\t\t|/;' $projectDirectory"matlab.SNP_analysis.log" >> $logName;
rm $outputName;
rm $projectDirectory"matlab.SNP_analysis.log"


##==============================================================================
## Generate final figures for dataset.
##------------------------------------------------------------------------------

echo "#==================================#" >> $logName;
echo "# Generate final combined figures. #" >> $logName;
echo "#==================================#" >> $logName;

echo "Generating final figures." >> $condensedLog;

echo "\t\tGenerating MATLAB script to generate combined CNV and SNP analysis figures from previous calculations." >> $logName;
outputName=$projectDirectory"processing4.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing4()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.final_figs.log');" >> $outputName;
echo "\tcd "$main_dir"Matlab/WGseq;" >> $outputName;
echo "\tanalyze_CNV_SNPs('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
echo "================================================================================================";
echo "== Final figures ===============================================================================";
echo "================================================================================================";
matlab -nosplash -r "run "$outputName"; exit;";
#chmod 0755 $outputName;
#chmod 0755 $projectDirectory"matlab.final_figs.log";
sed 's/^/\t\t|/;' $projectDirectory"matlab.final_figs.log" >> $logName;


##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------

echo "#=================================#" >> $logName;
echo "# Cleaning up intermediate files. #" >> $logName;
echo "#=================================#" >> $logName;

#rm $projectDirectory"";

## Generate "complete.txt" to indicate processing has completed normally.
completeFile=$projectDirectory"complete.txt";
echo "complete" > $completeFile;
echo "\tGenerated 'complete.txt' file." >> $logName;
