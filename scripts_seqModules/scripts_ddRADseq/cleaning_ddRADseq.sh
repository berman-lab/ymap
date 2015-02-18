#!/bin/bash -e
#
# cleaning_ddRADseq.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
main_dir=$(pwd)"/../../";


##==============================================================================
## Cleanup intermediate processing files.
##------------------------------------------------------------------------------
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";

echo "#==========================================#" >> $logName;
echo "# Cleaning up intermediate ddRADseq files. #" >> $logName;
echo "#==========================================#" >> $logName;
echo "Cleaning and archiving." >> $condensedLog;

if [ -f $projectDirectory"zipTemp.txt" ]
then
    rm $projectDirectory"zipTemp.txt";
	echo "\tzipTemp.txt" >> $logName;
fi
if [ -f $projectDirectory"processing1.m" ]
then
	rm $projectDirectory"processing1.m";
	echo "\tprocessing1.m" >> $logName;
fi
if [ -f $projectDirectory"processing2.m" ]
then
	rm $projectDirectory"processing2.m";
	echo "\tprocessing2.m" >> $logName;
fi
if [ -f $projectDirectory"processing3.m" ]
then
	rm $projectDirectory"processing3.m";
	echo "\tprocessing3.m" >> $logName;
fi
if [ -f $projectDirectory"processing4.m" ]
then
	rm $projectDirectory"processing4.m";
	echo "\tprocessing4.m" >> $logName;
fi
if [ -f $projectDirectory"matlab.ChARM.log" ]
then
	rm $projectDirectory"matlab.ChARM.log";
	echo "\tmatlab.ChARM.log" >> $logName;
fi
if [ -f $projectDirectory"matlab.CNV_and_GCbias.log" ]
then
	rm $projectDirectory"matlab.CNV_and_GCbias.log";
	echo "\tmatlab.CNV_and_GCbias.log" >> $logName;
fi
if [ -f $projectDirectory"matlab.final_figs.log" ]
then
	rm $projectDirectory"matlab.final_figs.log";
	echo "\tmatlat.final_figs.log" >> $logName;
fi
if [ -f $projectDirectory"matlab.SNP_analysis.log" ]
then
	rm $projectDirectory"matlab.SNP_analysis.log";
	echo "\tmatlab.SNP_analysis.log" >> $logName;
fi

if [ -f $projectDirectory"fig.ChARM_test.5.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.5.eps";
	echo "\tfig.ChARM_test.5.eps" >> $logName;
fi
if [ -f $projectDirectory"fig.ChARM_test.4.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.4.eps";
	echo "\tfig.ChARM_test.4.eps" >> $logName;
fi
if [ -f $projectDirectory"fig.ChARM_test.3.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.3.eps";
	echo "\tfig.ChARM_test.3.eps" >> $logName;
fi
if [ -f $projectDirectory"fig.ChARM_test.2.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.2.eps";
	echo "\tfig.ChARM_test.2.eps" >> $logName;
fi
if [ -f $projectDirectory"fig.ChARM_test.1.eps" ]
then
#	rm $projectDirectory"fig.ChARM_test.1.eps";
	echo "\tfig.ChARM_test.1.eps" >> $logName;
fi
if [ -f $projectDirectory"data_sorted.bam.bai" ]
then
#	rm $projectDirectory"data_sorted.bam.bai";
	echo "\tdata_sorted.bam.bai" >> $logName;
fi
if [ -f $projectDirectory"data_sorted.bam" ]
then
#	rm $projectDirectory"data_sorted.bam";
	echo "\tdata_sorted.bam" >> $logName;
fi
if [ -f $projectDirectory"data.pileup" ]
then
#	rm $projectDirectory"data.pileup";
	echo "\tdata.pileup" >> $logName;
fi
if [ -f $projectDirectory"data_indelRealigned.bam" ]
then
#	rm $projectDirectory"data_indelRealigned.bam";
	echo "\tdata_indelRealigned.bam" >> $logName;
fi
if [ -f $projectDirectory"data_indelRealigned.bai" ]
then
#	rm $projectDirectory"data_indelRealigned.bai";
	echo "\tdata_indelRealigned.bai" >> $logName;
fi
if [ -f $projectDirectory"data_forIndelRealigner.intervals" ]
then
#	rm $projectDirectory"data_forIndelRealigner.intervals";
	echo "\tdata_forIndelRealigner.intervals" >> $logName;
fi
if [ -f $projectDirectory"corrected_CNV.project.mat" ]
then
#	rm $projectDirectory"corrected_CNV.project.mat";
	echo "\tcorrected_CNV.project.mat" >> $logName;
fi
#	if [ -f $projectDirectory"fragment_CNV_data.mat" ]
#	then
#		rm $projectDirectory"fragment_CNV_data.mat";
#		echo "\tfragment_CNV_data.mat" >> $logName;
#	fi
if [ -f $projectDirectory"Common_ChARM.mat" ]
then
#	rm $projectDirectory"Common_ChARM.mat";
	echo "\tCommon_ChARM.mat" >> $logName;
fi
if [ -f $projectDirectory"Common_CNV.mat" ]
then
#	rm $projectDirectory"Common_CNV.mat";
	echo "\tCommon_CNV.mat" >> $logName;
fi
if [ -f $projectDirectory"SNP_v4.mat" ]
then
#	rm $projectDirectory"SNP_v4.mat";
	echo "\tSNP_v4.mat" >> $logName;
fi
#	if [ -f $projectDirectory"SNP_v4.all.mat" ]
#	then
#		rm $projectDirectory"SNP_v4.all.mat";
#		echo "\tSNP_v4.all.mat" >> $logName;
#	fi
if [ -f $projectDirectory"preprocessed_CNVs.ddRADseq.txt" ]
then
#	rm $projectDirectory"preprocessed_CNVs.ddRADseq.txt";
	echo "\tpreprocessed_CNVs.ddRADseq.txt" >> $logName;
fi
if [ -f $projectDirectory"preprocessed_SNPs.ddRADseq.CNV_filtered.txt" ]
then
#	rm $projectDirectory"preprocessed_SNPs.ddRADseq.CNV_filtered.txt";
	echo "\tpreprocessed_SNPs.ddRADseq.CNV_filtered.txt" >> $logName;
fi
if [ -f $projectDirectory"preprocessed_SNPs.ddRADseq.txt" ]
then
#	rm $projectDirectory"preprocessed_SNPs.ddRADseq.txt";
	echo "\tpreprocessed_SNPs.ddRADseq.txt" >> $logName;
fi
if [ -d $projectDirectory"fastqc_temp/" ]
then
	rm -rf $projectDirectory"fastqc_temp/";
	echo "\tfastqc_temp/" >> $logName;
fi
if [ -f $projectDirectory"allelic_ratios.txt" ]
then
#	rm $projectDirectory"allelic_ratios.txt";
	echo "\tallelic_ratios" >> $logName;
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
		if [ "$datafile1" != "null1" ]
		then
#			rm $projectDirectory$datafile1;
			echo "\t"$datafile1 >> $logName;
		fi
	else
		if [ "$datafile1" != "null1" ]
		then
#			rm $projectDirectory$datafile1;
			echo "\t"$datafile1 >> $logName;
		fi
		if [ "$datafile2" != "null2" ]
		then
#			rm $projectDirectory$datafile2;
			echo "\t"$datafile2 >> $logName;
		fi
	fi
#	rm $projectDirectory"datafiles.txt";
	echo "\tdatafiles.txt" >> $logName;
fi
if [ -f $projectDirectory"SNPdata_parent.txt" ]
then
#	rm $projectDirectory"SNPdata_parent.txt";
	echo "\tSNPdata_parent.txt" >> $logName;
fi


# Remove potential leftovers from upload restarts.
#pattern="*.zip";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.zip; fi
#echo "\t*.zip" >> $logName;
#pattern="*.gz";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.gz; fi
#echo "\t*.gz" >> $logName;
#pattern="*.bam";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.bam; fi
#echo "\t*.bam" >> $logName;
#pattern="*.sam";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.sam; fi
#echo "\t*.sam" >> $logName;
#pattern="*.fastq";
#if [ "$(echo $pattern)" != "$pattern" ]; then rm *.fastq; fi
#echo "\t*.fastq" >> $logName;

# Compress 'putative_SNPs_v1.txt' and 'SNP_CNVs_v1.txt'.
if [ -f $projectDirectory"putative_SNPs_v4.txt" ]
then
	zip -9 $projectDirectory"putative_SNPs_v4.zip" $projectDirectory"putative_SNPs_v4.txt";
#	rm $projectDirectory"putative_SNPs_v4.txt";
	echo "\tputative_SNPs_v4.txt" >> $logName;
fi
if [ -f $projectDirectory"SNP_CNV_v1.txt" ]
then
	zip -9 $projectDirectory"SNP_CNV_v1.zip" $projectDirectory"SNP_CNV_v1.txt";
#	rm $projectDirectory"SNP_CNV_v1.txt";
	echo "\tSNP_CNV_v1.txt" >> $logName;
fi

## Generate "complete.txt" to indicate processing has completed normally.
timestamp=$(date +%T);

	timesLogFile=$main_dir"completion_times.log";
	if [ -f $timesLogFile ]
	then
		echo -n $user"("$project")[WGseq " >> $timesLogFile;
		cat $projectDirectory"dataType.txt" >> $timesLogFile;
		echo -n "]\t" >> $timesLogFile;
		cat $projectDirectory"working.txt" >> $timesLogFile;
		echo " -> "$timestamp >> $timesLogFile;
	fi

completeFile=$projectDirectory"complete.txt";
echo "complete" > $completeFile;
echo "\tGenerated 'complete.txt' file." >> $logName;
chmod 0755 $completeFile;

if [ -f $projectDirectory"working.txt" ]
then
	mv $projectDirectory"working.txt" $projectDirectory"working_done.txt";
	echo "\tworking.txt" >> $logName;
fi

#rm $projectDirectory"process_log.txt";
#rm $projectDirectory"condensed_log.txt";
