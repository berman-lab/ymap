#!/bin/bash -e
#
# cleaning_WGseq.sh
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

. $main_dir"config.sh";
if [ $debug -eq 1 ];
then
	echo "\tReached cleanup stage, but skipping it because the debug flag is on." >> $logName;
	echo "\tCreating complete.txt, so that the front-end recognizes the completion." >> $logName;
	
	completeFile=$projectDirectory"complete.txt";
	echo "complete" > $completeFile;
	timestamp=$(date +%T);
	echo $timestamp >> $completeFile;
	echo "\tGenerated 'complete.txt' file." >> $logName;
	chmod 0744 $completeFile;
	
	## changing working.txt to working_done.txt
	if [ -f $projectDirectory"working.txt" ]
	then
		mv $projectDirectory"working.txt" $projectDirectory"working_done.txt";
		echo "\t changed working.txt to working_done.txt" >> $logName;
	fi
	
	exit 0;
fi

echo "#=======================================#" >> $logName;
echo "# Cleaning up intermediate WGseq files. #" >> $logName;
echo "#=======================================#" >> $logName;
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

if [ -f $projectDirectory"matlab.CNV_and_GCbias.log" ]
then
	rm $projectDirectory"matlab.CNV_and_GCbias.log";
	echo "\tmatlab.CNV_and_GCbias.log" >> $logName;
fi

if [ -f $projectDirectory"matlab.ChARM.log" ]
then
	rm $projectDirectory"matlab.ChARM.log";
	echo "\tmatlab.ChARM.log" >> $logName;
fi

if [ -f $projectDirectory"matlab.SNP_analysis.log" ]
then
	rm $projectDirectory"matlab.SNP_analysis.log";
	echo "\tmatlab.SNP_analysis.log" >> $logName;
fi

if [ -f $projectDirectory"matlab.final_figs.log" ]
then
	rm $projectDirectory"matlab.final_figs.log";
	echo "\tmatlab.final_figs.log" >> $logName;
fi

if [ -f $projectDirectory"data_sorted.bam.bai" ]
then
	rm $projectDirectory"data_sorted.bam.bai";
	echo "\tdata_sorted.bam.bai" >> $logName;
fi

if [ -f $projectDirectory"data_sorted.bam" ]
then
	rm $projectDirectory"data_sorted.bam";
	echo "\tdata_sorted.bam" >> $logName;
fi
if [ -f $projectDirectory"data.bam" ]
then
	rm $projectDirectory"data.bam";
	echo "\tdata.bam" >> $logName;
fi
if [ -f $projectDirectory"data.pileup" ]
then
	rm $projectDirectory"data.pileup";
	echo "\tdata.pileup" >> $logName;
fi
if [ -f $projectDirectory"data_indelRealigned.bam" ]
then
	rm $projectDirectory"data_indelRealigned.bam";
	echo "\tdata_indelRealigned.bam" >> $logName;
fi
if [ -f $projectDirectory"data_indelRealigned.bai" ]
then
	rm $projectDirectory"data_indelRealigned.bai";
	echo "\tdata_indelRealigned.bai" >> $logName;
fi

if [ -d $projectDirectory"fastqc_temp/" ]
then
	rm -rf $projectDirectory"fastqc_temp/";
	echo "\tfastqc_temp/" >> $logName;
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
		# deleting only if a valid file name is written
		if [ "$datafile1" != "null1" ]
		then
			if [ -f $projectDirectory$datafile1 ]
			then
				rm $projectDirectory$datafile1;
				echo "\t"$datafile1 >> $logName;
			fi
		fi
	else
		if [ -f $projectDirectory$datafile1 ]
		then
			rm $projectDirectory$datafile1;
			echo "\t"$datafile1 >> $logName;
		fi
		if [ -f $projectDirectory$datafile2 ]
		then
			rm $projectDirectory$datafile2;
			echo "\t"$datafile2 >> $logName;
		fi
	fi
	rm $projectDirectory"datafiles.txt";
	echo "\tdatafiles.txt" >> $logName;
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
	rm $projectDirectory"putative_SNPs_v4.txt";
	echo "\tputative_SNPs_v4.txt (created putative_SNPs_v4.zip)" >> $logName;
fi
if [ -f $projectDirectory"SNP_CNV_v1.txt" ]
then
	zip -9 $projectDirectory"SNP_CNV_v1.zip" $projectDirectory"SNP_CNV_v1.txt";
	rm $projectDirectory"SNP_CNV_v1.txt";
	echo "\tSNP_CNV_v1.txt (created SNP_CNV_v1.zip)" >> $logName;
fi

## Generate "complete.txt" to indicate processing has completed normally.
timestamp=$(date +%T);

	timesLogFile=$main_dir"completion_times.log";
	if [ -f $timesLogFile ]
	then
		echo -n $user"("$project")[WGseq " >> $timesLogFile;
		cat $projectDirectory"dataFormat.txt" >> $timesLogFile;
		echo -n "]\t" >> $timesLogFile;
		cat $projectDirectory"working.txt" >> $timesLogFile;
		echo " -> "$timestamp >> $timesLogFile;
	fi

completeFile=$projectDirectory"complete.txt";
echo "complete" > $completeFile;
echo $timestamp >> $completeFile;
echo "\tGenerated 'complete.txt' file." >> $logName;
chmod 0755 $completeFile;

if [ -f $projectDirectory"working.txt" ]
then
	mv $projectDirectory"working.txt" $projectDirectory"working_done.txt";
	echo "\tworking.txt" >> $logName;
fi
