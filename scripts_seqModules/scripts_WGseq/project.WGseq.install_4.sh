#!/bin/bash -e
#
# project.WGseq.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
main_dir=$(pwd)"/../../";

echo "";
echo "Input to : project.WGseq.install_4.sh";
echo "\tuser     = "$user;
echo "\tproject  = "$project;
echo "\tmain_dir = "$main_dir;
echo "";

# load local installed program location variables.
. $main_dir"local_installed_programs.sh";

##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------

projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";

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

# Get parent name from "parent.txt" in project directory.
projectParent=$(head -n 1 $projectDirectory"parent.txt");
echo "\tparentProject = '"$projectParent"'" >> $logName;

# Determine location of project being used.
if [ -d $main_dir"users/"$user"/projects/"$projectParent"/" ]
then
    projectParentDirectory=$main_dir"users/"$user"/projects/"$projectParent"/";
    projectParentUser=$user;
elif [ -d $main_dir"users/default/projects/"$projectParent"/" ]
then
    projectParentDirectory=$main_dir"users/default/projects/"$projectParent"/";
    projectParentUser="default";
fi
echo "\tprojectParentDirectory = '"$projectParentDirectory"'" >> $logName;



##==============================================================================
## Perform CGH analysis, with GC-correction, on dataset.
##------------------------------------------------------------------------------
echo "#==========================#" >> $logName;
echo "# CGH analysis of dataset. #" >> $logName;
echo "#==========================#" >> $logName;
echo "Preprocessing CNV data.   (~10 min for 1.6 Gbase genome dataset.)" >> $condensedLog;

if [ -f $projectDirectory"preprocessed_CNVs.txt" ]
then
	echo "\t\tCNV data already preprocessed with python script : 'scripts_seqModules/scripts_WGseq/dataset_process_for_CNV_analysis.WGseq.py'" >> $logName;
else
	echo "\t\tPreprocessing CNV data with python script : 'scripts_seqModules/scripts_WGseq/dataset_process_for_CNV_analysis.WGseq.py'" >> $logName;
	$python_exec $main_dir"scripts_seqModules/scripts_WGseq/dataset_process_for_CNV_analysis.WGseq.py" $user $project $genome $genomeUser $main_dir $logName  > $projectDirectory"preprocessed_CNVs.txt";
	echo "\t\tpre-processing complete." >> $logName;
fi

echo "Analyzing and mapping CNVs." >> $condensedLog;

echo "\t\tGenerating MATLAB script to perform CNV analysis of dataset, with GC-correction." >> $logName;
outputName=$projectDirectory"processing1.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing1()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.CNV_and_GCbias.log');" >> $outputName;
echo "\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $outputName;
echo "\tanalyze_CNVs_1('$main_dir','$user','$genomeUser','$project','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t|\tfunction [] = processing1()" >> $logName;
echo "\t|\t\tdiary('"$projectDirectory"matlab.CNV_and_GCbias.log');" >> $logName;
echo "\t|\t\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $logName;
echo "\t|\t\tanalyze_CNVs_1('$main_dir','$user','$genomeUser','$project','$genome','$ploidyEstimate','$ploidyBase');" >> $logName;
echo "\t|\tend" >> $logName;

echo "\t\tCalling MATLAB." >> $logName;
$matlab_exec -nosplash -r "run "$outputName"; exit;";
echo "\t\tMATLAB log from CNV analysis." >> $logName;
sed 's/^/\t\t\t|/;' $projectDirectory"matlab.CNV_and_GCbias.log" >> $logName;


##==============================================================================
## Perform ChARM analysis of dataset.
##------------------------------------------------------------------------------
echo "#============================#" >> $logName;
echo "# ChARM analysis of dataset. #" >> $logName;
echo "#============================#" >> $logName;
echo "Analyzing CNV edges." >> $condensedLog;

if [ -f $projectDirectory"Common_ChARM.mat" ]
then
	echo "\t\tChARM analysis already completed." >> $logName;
else
	echo "\t\tGenerating MATLAB script to perform ChARM analysis of dataset." >> $logName;
	outputName=$projectDirectory"processing2.m";
	echo "\t\toutputName = "$outputName >> $logName;

	echo "function [] = processing2()" > $outputName;
	echo "\tdiary('"$projectDirectory"matlab.ChARM.log');" >> $outputName;
	echo "\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $outputName;
	echo "\tChARM_v4('$project','$user','$genome','$genomeUser','$main_dir');" >> $outputName;
	echo "end" >> $outputName;

	echo "\t|\tfunction [] = processing2()" >> $logName;
	echo "\t|\t\tdiary('"$projectDirectory"matlab.ChARM.log');" >> $logName;
	echo "\t|\t\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $logName;
	echo "\t|\t\tChARM_v4('$project','$user','$genome','$genomeUser','$main_dir');" >> $logName;
	echo "\t|\tend" >> $logName;

	echo "\t\tCalling MATLAB." >> $logName;
	echo "================================================================================================";
	echo "== ChARM analysis ==============================================================================";
	echo "================================================================================================";
	$matlab_exec -nosplash -r "run "$outputName"; exit;";
	echo "\t\tMATLAB log from ChARM analysis." >> $logName;
	sed 's/^/\t\t\t|/;' $projectDirectory"matlab.ChARM.log" >> $logName;
fi

##==============================================================================
## Perform SNP/LOH analysis on dataset.
##------------------------------------------------------------------------------
if [ "$project" = "$projectParent" ]
then
	echo "#==========================#" >> $logName;
	echo "# SNP analysis of dataset. #" >> $logName;
	echo "#==========================#" >> $logName;
	echo "Preprocessing SNP data.   (~20 min for SNP analysis of 1.6 Gbase genome dataset.)" >> $condensedLog;
else
	echo "#==========================#" >> $logName;
	echo "# LOH analysis of dataset. #" >> $logName;
	echo "#==========================#" >> $logName;
	echo "Preprocessing SNP data.   (~4 hrs for LOH analysis of 1.6 Gbase genome dataset.)" >> $condensedLog;
fi

if [ -f $projectDirectory"preprocessed_SNPs.txt" ]
then
	echo "\t\tSNP data already preprocessed with python script : 'scripts_seqModules/scripts_WGseq/dataset_process_for_SNP_analysis.WGseq.py'" >> $logName;
else
	echo "\t\tPreprocessing SNP data with python script : 'scripts_seqModules/scripts_WGseq/dataset_process_for_SNP_analysis.WGseq.py'" >> $logName;
	if [ -f $projectParentDirectory"putative_SNPs_v4.txt" ]
	then
		echo "\t\tParent SNP data already decompressed." >> $logName;
		cp $projectParentDirectory"putative_SNPs_v4.txt" $projectDirectory"SNPdata_parent.txt";
	else
		echo "\t\tDecompressing parent SNP data." >> $logName;
		unzip -j $projectParentDirectory"putative_SNPs_v4.zip" -d $projectDirectory;
        mv $projectDirectory"putative_SNPs_v4.txt" $projectDirectory"SNPdata_parent.txt";
	fi

	# preprocess parent for comparison. abbey
	$python_exec $main_dir"scripts_seqModules/scripts_hapmaps/hapmap.preprocess_parent.py" $genome $genomeUser $project $user $projectParent $projectParentUser $main_dir LOH > $projectDirectory"SNPdata_parent.temp.txt";

	rm $projectDirectory"SNPdata_parent.txt";
	mv $projectDirectory"SNPdata_parent.temp.txt" $projectDirectory"SNPdata_parent.txt";

	$python_exec $main_dir"scripts_seqModules/scripts_WGseq/dataset_process_for_SNP_analysis.WGseq.py" $genome $genomeUser $projectParent $projectParentUser $project $user $main_dir $logName LOH > $projectDirectory"preprocessed_SNPs.txt";
	echo "\t\tpre-processing complete." >> $logName;
fi

echo "Mapping SNPs." >> $condensedLog;
echo "\t\tGenerating MATLAB script to perform SNP analysis of dataset." >> $logName;
outputName=$projectDirectory"processing3.m";
echo "\t\toutputName = "$outputName >> $logName;

echo "function [] = processing3()" > $outputName;
echo "\tdiary('"$projectDirectory"matlab.SNP_analysis.log');" >> $outputName;
echo "\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $outputName;
echo "\tanalyze_SNPs_hapmap('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t|\tfunction [] = processing3()" >> $logName;
echo "\t|\t\tdiary('"$projectDirectory"matlab.SNP_analysis.log');" >> $logName;
echo "\t|\t\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $logName;
echo "\t|\t\tanalyze_SNPs_hapmap('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $logName;
echo "\t|\tend" >> $logName;

echo "\t\tCalling MATLAB." >> $logName;
echo "================================================================================================";
echo "== SNP analysis ================================================================================";
echo "================================================================================================";
$matlab_exec -nosplash -r "run "$outputName"; exit;";
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
echo "\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $outputName;
echo "\tanalyze_CNV_SNPs_hapmap('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $outputName;
echo "end" >> $outputName;

echo "\t|\tfunction [] = processing4()" >> $logName;
echo "\t|\t\tdiary('"$projectDirectory"matlab.final_figs.log');" >> $logName;
echo "\t|\t\tcd "$main_dir"scripts_seqModules/scripts_WGseq;" >> $logName;
echo "\t|\t\tanalyze_CNV_SNPs_hapmap('$main_dir','$user','$genomeUser','$project','$projectParent','$genome','$ploidyEstimate','$ploidyBase');" >> $logName;
echo "\t|\tend" >> $logName;

echo "\t\tCalling MATLAB.   (Log will be appended here after completion.)" >> $logName;
echo "================================================================================================";
echo "== Final figures ===============================================================================";
echo "================================================================================================";
$matlab_exec -nosplash -r "run "$outputName"; exit;";
echo "\t\tMATLAB log from final figure generation." >> $logName;
sed 's/^/\t\t|/;' $projectDirectory"matlab.final_figs.log" >> $logName;
echo "finished all processing, moving to Cleaning up intermediate WGseq files" >> $condensedLog;


##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------
echo "running: " $main_dir"scripts_seqModules/scripts_WGseq/cleaning_WGseq.sh" $user $project $main_dir >> $logName;
sh $main_dir"scripts_seqModules/scripts_WGseq/cleaning_WGseq.sh" $user $project $main_dir 2>> $logName;
