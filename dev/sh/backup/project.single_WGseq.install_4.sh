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
user=$1;
project=$2;
main_dir=$3;


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
matlab -nosplash -r "run "$outputName";";
chmod 0755 $outputName;
chmod 0755 $projectDirectory"matlab.CNV_and_GCbias.log";
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
echo "Analyzing CNV edges." >> $condensedLog;

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
matlab -nosplash -r "run "$outputName";";
chmod 0755 $outputName;
chmod 0755 $projectDirectory"matlab.ChARM.log";
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
matlab -nosplash -r "run "$outputName";";
chmod 0755 $outputName;
chmod 0755 $projectDirectory"matlab.SNP_analysis.log";
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
echo "== Final Figures ===============================================================================";
echo "================================================================================================";
matlab -nosplash -r "run "$outputName";";
chmod 0755 $outputName;
chmod 0755 $projectDirectory"matlab.final_figs.log";
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
