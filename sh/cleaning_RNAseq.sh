#!/bin/bash -e
#
# cleaning_RNAseq.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
main_dir=$3;

#user="darren";
#project="BG2_test_FF";
#main_dir="/heap/hapmap/bermanlab/";

##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------

projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";

echo "#========================================#" >> $logName;
echo "# Cleaning up intermediate RNAseq files. #" >> $logName;
echo "#========================================#" >> $logName;

if [ -f $projectDirectory"zipTemp.txt" ]
then
	rm $projectDirectory"zipTemp.txt";
fi
if [ -f $projectDirectory"txt.CNV-map.3.txt" ]
then
#	rm $projectDirectory"txt.CNV-map.3.txt";
fi
if [ -f $projectDirectory"showAnnotations.txt" ]
then
#	rm $projectDirectory"showAnnotations.txt";
fi
if [ -f $projectDirectory"processing1.m" ]
then
	rm $projectDirectory"processing1.m";
fi
if [ -f $projectDirectory"processing2.m" ]
then
	rm $projectDirectory"processing2.m";
fi
if [ -f $projectDirectory"processing3.m" ]
then
	rm $projectDirectory"processing3.m";
fi
if [ -f $projectDirectory"processing4.m" ]
then
	rm $projectDirectory"processing4.m";
fi
if [ -f $projectDirectory"ploidy.txt" ]
then
#	rm $projectDirectory"ploidy.txt";
fi
if [ -f $projectDirectory"matlab.CNV_and_GCbias.log" ]
then
	rm $projectDirectory"matlab.CNV_and_GCbias.log";
fi
if [ -f $projectDirectory"matlab.ChARM.log" ]
then
	rm $projectDirectory"matlab.ChARM.log";
fi
if [ -f $projectDirectory"matlab.SNP_analysis.log" ]
then
	rm $projectDirectory"matlab.SNP_analysis.log";
fi
if [ -f $projectDirectory"matlab.final_figs.log" ]
then
	rm $projectDirectory"matlab.final_figs.log";
fi
if [ -f $projectDirectory"fig.ChARM_test.5.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.5.eps";
fi
if [ -f $projectDirectory"fig.ChARM_test.4.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.4.eps";
fi
#if [ -f $projectDirectory"fig.ChARM_test.4.png" ]
#then
#	rm $projectDirectory"fig.ChARM_test.4.png";
#fi
if [ -f $projectDirectory"fig.ChARM_test.3.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.3.eps";
fi
#if [ -f $projectDirectory"fig.ChARM_test.3.png" ]
#then
#	rm $projectDirectory"fig.ChARM_test.3.png";
#fi
if [ -f $projectDirectory"fig.ChARM_test.2.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.2.eps";
fi
#if [ -f $projectDirectory"fig.ChARM_test.2.png" ]
#then
#	rm $projectDirectory"fig.ChARM_test.2.png";
#fi
if [ -f $projectDirectory"fig.ChARM_test.1.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.1.eps";
fi
#if [ -f $projectDirectory"fig.ChARM_test.1.png" ]
#then
#	rm $projectDirectory"fig.ChARM_test.1.png";
#fi
if [ -f $projectDirectory"fig.GCratio_vs_SNP.png" ]
then
#	rm $projectDirectory"fig.GCratio_vs_SNP.png";
fi
if [ -f $projectDirectory"data_sorted.bam.bai" ]
then
#	rm $projectDirectory"data_sorted.bam.bai";
fi
if [ -f $projectDirectory"data_sorted.bam" ]
then
#	rm $projectDirectory"data_sorted.bam";
fi
if [ -f $projectDirectory"data.pileup" ]
then
#	rm $projectDirectory"data.pileup";
fi
if [ -f $projectDirectory"data_indelRealigned.bam" ]
then
#	rm $projectDirectory"data_indelRealigned.bam";
fi
if [ -f $projectDirectory"data_indelRealigned.bai" ]
then
#	rm $projectDirectory"data_indelRealigned.bai";
fi
if [ -f $projectDirectory"data_forIndelRealigner.intervals" ]
then
#	rm $projectDirectory"data_forIndelRealigner.intervals";
fi
if [ -f $projectDirectory"data.bam" ]
then
#	rm $projectDirectory"data.bam";
fi
if [ -f $projectDirectory"CNV_v1.mat" ]
then
#	rm $projectDirectory"CNV_v1.mat";
fi
if [ -f $projectDirectory"Common_ChARM.mat" ]
then
#	rm $projectDirectory"Common_ChARM.mat";
fi
if [ -f $projectDirectory"Common_CNV.mat" ]
then
#	rm $projectDirectory"Common_CNV.mat";
fi
if [ -f $projectDirectory"SNP_v4.mat" ]
then
#	rm $projectDirectory"SNP_v4.mat";
fi
if [ -f $projectDirectory"preprocessed_CNVs.txt" ]
then
#	rm $projectDirectory"preprocessed_CNVs.txt";
fi
if [ -f $projectDirectory"preprocessed_SNPs.txt" ]
then
#	rm $projectDirectory"preprocessed_SNPs.txt";
fi
if [ -d $projectDirectory"fastqc_temp/" ]
then
#	rm -rf $projectDirectory"fastqc_temp/";
fi
if [ -f $projectDirectory"datafiles.txt" ]
then
	# Get first data file name from "datafiles.txt";
	datafile1=$(head -n 1 $projectDirectory"datafiles.txt");
	echo "\tdatafile 1 = '"$datafile1"'" >> $logName;
	# Get second data file name from "datafiles.txt";
	datafile2=$(tail -n 1 $projectDirectory"datafiles.txt");
	echo "\tdatafile 2 = '"$datafile2"'" >> $logName;
	if [ "$datafile1" = "$datafile2" ]
	then
		if [ "$datafile1" = "null1" ]
		then
#			rm $projectDirectory$datafile1;
		fi
	else
#		rm $projectDirectory$datafile1;
#		rm $projectDirectory$datafile2;
	fi
#	rm $projectDirectory"datafiles.txt";
fi
if [ -f $projectDirectory"SNPdata_parent.txt" ]
then
#	rm $projectDirectory"SNPdata_parent.txt";
fi


# Remove potential leftovers from upload restarts.
#pattern="*.zip";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.zip; fi
#pattern="*.gz";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.gz; fi
#pattern="*.bam";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.bam; fi
#pattern="*.sam";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.sam; fi
#pattern="*.fastq";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.fastq; fi


# Compress 'putative_SNPs_v1.txt' and 'SNP_CNVs_v1.txt'.
if [ -f $projectDirectory"putative_SNPs_v4.txt" ]
then
	zip -9 $projectDirectory"putative_SNPs_v4.zip" $projectDirectory"putative_SNPs_v4.txt";
#	rm $projectDirectory"putative_SNPs_v4.txt";
fi
if [ -f $projectDirectory"SNP_CNV_v1.txt" ]
then
	zip -9 $projectDirectory"SNP_CNV_v1.zip" $projectDirectory"SNP_CNV_v1.txt";
#	rm $projectDirectory"SNP_CNV_v1.txt";
fi

## Generate "complete.txt" to indicate processing has completed normally.
completeFile=$projectDirectory"complete.txt";
echo "complete" > $completeFile;
echo "\tGenerated 'complete.txt' file." >> $logName;
chmod 0755 $completeFile;

if [ -f $projectDirectory"working.txt" ]
then
	rm $projectDirectory"working.txt";
fi

#rm $projectDirectory"process_log.txt";
#rm $projectDirectory"condensed_log.txt";
