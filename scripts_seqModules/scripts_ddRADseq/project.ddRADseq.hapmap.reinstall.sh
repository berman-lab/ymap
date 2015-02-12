#!/bin/bash -e
#
# project.ddRADseq.hapmap.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
hapmap=$3;
main_dir=$(pwd)"/";

user='default';
hapmap='C_albicans_SC5314_A21-s02-m09-r07';

project='Fig_06C.YJB12712-d2';



##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";


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

# Determine location of hapmap being used.
echo "\thapmap = '$hapmap'" >> $logName;
hapmapTestDir1=$main_dir"users/"$user"/hapmaps/"$hapmap"/";
echo "looking for dir '$hapmapTestDir1'" >> $logName;
echo "Script executed from: ${PWD}" >> $logName;
testingDir=$main_dir"users/default/hapmaps/"$hapmap"/";
if [ -d $testing ]
then
	echo "testing dir '$testingDir' found." >> $logName;
else
	echo "testing dir '$testingDir' not found." >> $logName;
fi

testingDir1=$main_dir"users/"$user"/hapmaps/"$hapmap"/";
testingDir2=$main_dir"users/default/hapmaps/"$hapmap"/";
echo "testingDir1 = '$testingDir1'" >> $logName;
echo "testingDir2 = '$testingDir2'" >> $logName;
if [ -d $testingDir1 ]
then
    hapmapDirectory=$main_dir"users/"$user"/hapmaps/"$hapmap"/";
    hapmapUser=$user;
    hapmapUsed=1;
	echo "hapmap found in user dir." >> $logName;
elif [ -d $testingDir2 ]
then
    hapmapDirectory=$main_dir"users/default/hapmaps/"$hapmap"/";
    hapmapUser="default";
    hapmapUsed=1;
	echo "hapmap found in default dir." >> $logName;
else
    hapmapUsed=0;
	echo "hapmap not found." >> $logName;
fi
echo "\thapmapDirectory = '"$hapmapDirectory"'" >> $logName;

cp $hapmapDirectory"colors.txt" $projectDirectory"colors.txt";


reflocation=$main_dir"users/"$genomeUser"/genomes/"$genome"/";                 # Directory where FASTA file is kept.
FASTA=`sed -n 1,1'p' $reflocation"reference.txt"`;                             # Name of FASTA file.
FASTAname=$(echo $FASTA | sed 's/.fasta//g');                                  # name of genome file, without file type.
RestrctionEnzymes=`sed -n 1,1'p' $projectDirectory"restrictionEnzymes.txt"`;   # Name of restriction enxyme list file.
ddRADseq_FASTA=$FASTAname"."$RestrctionEnzymes".fasta";                        # Name of digested reference for ddRADseq analysis, using chosen restriction enzymes.


##==============================================================================
## Preprocess ddRADseq CNV information.
##------------------------------------------------------------------------------
if [ -f $projectDirectory"preprocessed_CNVs.ddRADseq.txt" ]
then
	echo "\t\tCNV data already preprocessed with python script : 'py/dataset_process_for_CNV_analysis.ddRADseq.py'" >> $logName;
else
	echo "\t\tPreprocessing CNV data with python script : 'py/dataset_process_for_CNV_analysis.ddRADseq.py'" >> $logName;
	python $main_dir"py/dataset_process_for_CNV_analysis.ddRADseq.py" $user $project $genome $genomeUser $main_dir $RestrctionEnzymes $logName  > $projectDirectory"preprocessed_CNVs.ddRADseq.txt";
	echo "\t\tpre-processing complete." >> $logName;
fi


##==============================================================================
## Preprocess ddRADseq SNP information.
##------------------------------------------------------------------------------
if [ -f $projectDirectory"preprocessed_SNPs.ddRADseq.txt" ]
then
	echo "\t\tSNP data already preprocessed with python script : 'py/dataset_process_for_SNP_analysis.ddRADseq.py'" >> $logName;
else
	echo "\t\tPreprocessing SNP data with python script : 'py/dataset_process_for_SNP_analysis.ddRADseq.py'" >> $logName;
	if [ -f $parentDirectory"putative_SNPs_v4.txt" ]
	then
		echo "\t\tParent SNP data already decompressed." >> $logName;
		cp $parentDirectory"putative_SNPs_v4.txt" $projectDirectory"SNPdata_parent.txt";
	else
		echo "\t\tDecompressing parent SNP data." >> $logName;
		cd $parentDirectory;
		unzip -j putative_SNPs_v4.zip;
		cd $main_dir;
		cp $parentDirectory"putative_SNPs_v4.txt" $projectDirectory"SNPdata_parent.txt";
	fi

	# preprocess parent for comparison. abbey
	python $main_dir"py/hapmap.preprocess_parent.py" $genome $genomeUser $project $user $parent $parentUser $main_dir LOH > $projectDirectory"SNPdata_parent.temp.txt";

	rm $projectDirectory"SNPdata_parent.txt";
	mv $projectDirectory"SNPdata_parent.temp.txt" $projectDirectory"SNPdata_parent.txt";

	python $main_dir"py/dataset_process_for_SNP_analysis.ddRADseq.py" $genome $genomeUser $parent $parentUser $project $user $main_dir $RestrctionEnzymes $logName LOH > $projectDirectory"preprocessed_SNPs.ddRADseq.txt";
	chmod 0777 $projectDirectory"preprocessed_SNPs.ddRADseq.txt";
	echo "\t\tpre-processing complete." >> $logName;
fi


##==============================================================================
## Perform CGH analysis, with GC-correction, on dataset.
##------------------------------------------------------------------------------
echo "#==========================#" >> $logName;
echo "# CGH analysis of dataset. #" >> $logName;
echo "#==========================#" >> $logName;
echo "Preprocessing CNV data.   (~10 min for 1.6 Gbase genome dataset.)" >> $condensedLog;

if [ -f $projectDirectory"corrected_CNV.project.mat" ]
then
	echo "\t\tCNV analysis already complete." >> $logName;
else
	echo "Analyzing and mapping CNVs." >> $condensedLog;

	echo "\t\tGenerating MATLAB script to perform CNV analysis of dataset, with GC-correction." >> $logName;
	outputName=$projectDirectory"processing1.m";
	echo "\t\toutputName = "$outputName >> $logName;

	echo "function [] = processing1()" > $outputName;
	echo "\tdiary('"$projectDirectory"matlab.CNV_and_GCbias.log');" >> $outputName;
	echo "\tcd "$main_dir"Matlab/ddRADseq;" >> $outputName;
	echo "\tanalyze_CNVs_RADseq_3('$main_dir','$user','$genomeUser','$project','$parent','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
	echo "end" >> $outputName;

	echo "\t|\tfunction [] = processing1()" >> $logName;
	echo "\t|\t\tdiary('"$projectDirectory"matlab.CNV_and_GCbias.log');" >> $logName;
	echo "\t|\t\tcd "$main_dir"Matlab/ddRADseq;" >> $logName;
	echo "\t|\t\tanalyze_CNVs_RADseq_3('$main_dir','$user','$genomeUser','$project','$parent','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $logName;
	echo "\t|\tend" >> $logName;

	echo "\t\tCalling MATLAB." >> $logName;
	matlab -nosplash -r "run "$outputName"; exit;";
	echo "\t\tMATLAB log from CNV analysis." >> $logName;
	sed 's/^/\t\t\t|/;' $projectDirectory"matlab.CNV_and_GCbias.log" >> $logName;
fi


##==============================================================================
## Perform ChARM analysis of dataset.
##------------------------------------------------------------------------------
echo "#============================#" >> $logName;
echo "# ChARM analysis of dataset. #" >> $logName;
echo "#============================#" >> $logName;
echo "Analyzing CNV edges." >> $condensedLog;

if [ -f $projectDirectory"Common_ChARM.mat" ]
then
	echo "\t\tChARM analysis already complete." >> $logName;
else
	echo "\t\tGenerating MATLAB script to perform ChARM analysis of dataset." >> $logName;
	outputName=$projectDirectory"processing2.m";
	echo "\t\toutputName = "$outputName >> $logName; 

	echo "function [] = processing2()" > $outputName;
	echo "\tdiary('"$projectDirectory"matlab.ChARM.log');" >> $outputName;
	echo "\tcd "$main_dir"Matlab/ChARM;" >> $outputName;
	echo "\tChARM_v4('$project','$user','$genome','$genomeUser','$main_dir');" >> $outputName;
	echo "end" >> $outputName;

	echo "\t|\tfunction [] = processing2()" >> $logName;
	echo "\t|\t\tdiary('"$projectDirectory"matlab.ChARM.log');" >> $logName;
	echo "\t|\t\tcd "$main_dir"Matlab/ChARM;" >> $logName;
	echo "\t|\t\tChARM_v4('$project','$user','$genome','$genomeUser','$main_dir');" >> $logName;
	echo "\t|\tend" >> $logName;

	echo "\t\tCalling MATLAB." >> $logName;
	echo "================================================================================================";
	echo "== ChARM analysis ==============================================================================";
	echo "================================================================================================";
	matlab -nosplash -r "run "$outputName"; exit;";
	echo "\t\tMATLAB log from ChARM analysis." >> $logName;
	sed 's/^/\t\t\t|/;' $projectDirectory"matlab.ChARM.log" >> $logName;
fi


##==============================================================================
## Perform SNP/LOH analysis on dataset.   ...must be redone for ddRADseq, specificially.
##------------------------------------------------------------------------------
if [ "$project" = "$parent" ]
then
    echo "#==========================#" >> $logName;
    echo "# LOH analysis of dataset. #" >> $logName;
    echo "#==========================#" >> $logName;
    echo "Preprocessing SNP data.   (~4 hrs for LOH analysis of 1.6 Gbase genome dataset.)" >> $condensedLog;
else
    echo "#==========================#" >> $logName;
    echo "# SNP analysis of dataset. #" >> $logName;
    echo "#==========================#" >> $logName;
    echo "Preprocessing SNP data.   (~20 min for SNP analysis of 1.6 Gbase genome dataset.)" >> $condensedLog;
fi

echo "Mapping SNPs." >> $condensedLog;
echo "\t\tGenerating MATLAB script to perform SNP analysis of dataset." >> $logName;
outputName=$projectDirectory"processing3.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing3()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.SNP_analysis.log');" >> $outputName;
echo "\tcd "$main_dir"Matlab/ddRADseq;" >> $outputName;
echo "\tanalyze_SNPs_RADseq('$main_dir','$user','$genomeUser','$project','$parent','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t|\tfunction [] = processing3()" >> $logName;
echo "\t|\t\tdiary('"$projectDirectory"matlab.SNP_analysis.log');" >> $logName;
echo "\t|\t\tcd "$main_dir"Matlab/ddRADseq;" >> $logName;
echo "\t|\t\tanalyze_SNPs_RADseq('$main_dir','$user','$genomeUser','$project','$parent','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $logName;
echo "\t|\tend" >> $logName;

echo "\t\tCalling MATLAB." >> $logName;
echo "================================================================================================";
echo "== SNP analysis ================================================================================";
echo "================================================================================================";
matlab -nosplash -r "run "$outputName"; exit;";
echo "\t\tMATLAB log from SNP analysis." >> $logName;
sed 's/^/\t\t\t|/;' $projectDirectory"matlab.SNP_analysis.log" >> $logName;


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
echo "\tcd "$main_dir"Matlab/ddRADseq;" >> $outputName;
echo "\tanalyze_CNV_SNPs_RADseq('$main_dir','$user','$genomeUser','$project','$parent','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t|\tfunction [] = processing4()" >> $logName;
echo "\t|\t\tdiary('"$projectDirectory"matlab.final_figs.log');" >> $logName;
echo "\t|\t\tcd "$main_dir"Matlab/ddRADseq;" >> $logName;
echo "\t|\t\tanalyze_CNV_SNPs_RADseq('$main_dir','$user','$genomeUser','$project','$parent','$hapmap','$genome','$ploidyEstimate','$ploidyBase');" >> $logName;
echo "\t|\tend" >> $logName;

echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
matlab -nosplash -r "run "$outputName"; exit;";
sed 's/^/\t\t|/;' $projectDirectory"matlab.final_figs.log" >> $logName;


##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------
sh $main_dir"sh/cleaning_ddRADseq.sh" $user $project $main_dir;
