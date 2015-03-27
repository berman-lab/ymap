#!/bin/bash -e
#
# project.WGseq.hapmap.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;


user='darren1';
project='SC5314_WGseq_';
hapmap='';

main_dir=$(pwd)"/../../";


##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";


if [ -z "$hapmap"]
then
	# Get parent name used from project's "parent.txt" file.
	parent=$(head -n 1 $projectDirectory"parent.txt");
	echo "\tparent = '"$parent"'" >> $logName;
	# Determine location of parent.
	if [ -d $main_dir"users/"$user"/projects/"$parent"/" ]
	then
		parentDirectory=$main_dir"users/"$user"/projects/"$parent"/";
		parentUser=$user;
	elif [ -d $main_dir"users/default/projects/"$parent"/" ]
	then
		parentDirectory=$main_dir"users/default/projects/"$parent"/";
		parentUser="default";
	fi
	echo "\tparentDirectory = '"$parentDirectory"'" >> $logName;
fi


# Get genome name used from project's "genome.txt" file.
genome=$(head -n 1 $projectDirectory"genome.txt");
echo "\tgenome = '"$genome"'" >> $logName;
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

# Get ploidy estimate from "ploidy.txt" in project directory.
ploidyEstimate=$(head -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyEstimate = '"$ploidyEstimate"'" >> $logName;

# Get ploidy baseline from "ploidy.txt" in project directory.
ploidyBase=$(tail -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyBase = '"$ploidyBase"'" >> $logName;

if [ -z "$hapmap" ]
then
	hapmapUsed=0;
else
	# Determine location of hapmap being used.
	if [ -d $main_dir"users/"$user"/hapmaps/"$hapmap"/" ]
	then
		hapmapDirectory=$main_dir"users/"$user"/hapmaps/"$hapmap"/";
		hapmapUser=$user;
		hapmapUsed=1
	elif [ -d $main_dir"users/default/hapmaps/"$hapmap"/" ]
	then
		hapmapDirectory=$main_dir"users/default/hapmaps/"$hapmap"/";
		hapmapUser="default";
		hapmapUsed=1;
	else
		hapmapUsed=0;
	fi
	echo "\thapmapDirectory = '"$hapmapDirectory"'" >> $logName;
	cp $hapmapDirectory"colors.txt" $projectDirectory"colors.txt";
fi


reflocation=$main_dir"users/"$genomeUser"/genomes/"$genome"/";                 # Directory where FASTA file is kept.
FASTA=`sed -n 1,1'p' $reflocation"reference.txt"`;                             # Name of FASTA file.
FASTAname=$(echo $FASTA | sed 's/.fasta//g');                                  # name of genome file, without file type.


##==============================================================================
## Generate script to re-run terminal visualization Matlab code.
##------------------------------------------------------------------------------
echo "#======================================#" >> $logName;
echo "# Re-perform visualization of dataset. #" >> $logName;
echo "#======================================#" >> $logName;

echo "\tGenerating MATLAB script to perform CNV analysis of dataset, with GC-correction." >> $logName;
outputName=$projectDirectory"processing_Rerun.m";
echo "\toutputName = "$outputName >> $logName;

echo "function [] = processing_Rerun()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.rerun_visualization.log');" >> $outputName;
echo "\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $outputName;

echo     "\tanalyze_CNVs_1(         '$main_dir','$user','$genomeUser','$project',          '$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
if [ -z "$hapmap" ]
then
	echo "\tanalyze_SNPs_hapmap(    '$main_dir','$user','$genomeUser','$project','$parent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
	echo "\tanalyze_CNV_SNPs_hapmap('$main_dir','$user','$genomeUser','$project','$parent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
else
	echo "\tanalyze_SNPs_hapmap(    '$main_dir','$user','$genomeUser','$project','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
	echo "\tanalyze_CNV_SNPs_hapmap('$main_dir','$user','$genomeUser','$project','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
fi

echo "end" >> $outputName;

echo "\tCalling MATLAB." >> $logName;
matlab -nosplash -nodesktop -r "run "$outputName"; exit;";
echo "\tMATLAB log from redo of visualization.." >> $logName;
sed 's/^/\t\t|/;' $projectDirectory"matlab.rerun_visualization.log" >> $logName;
