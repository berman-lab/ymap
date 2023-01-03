#!/bin/bash -e
#
# project.single_IonExpressSeq.install_3.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
main_dir=$(pwd)"/";

#user="darren";
#project="morlurie_30-26_b";
#main_dir="/heap/hapmap/bermanlab/";

##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------

# Setup process_log.txt file.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";
echo "#.............................................................................." >> $logName;

# import locations of auxillary software for pipeline analysis.
. $main_dir"local_installed_programs.sh";

# Define project directory.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";

echo "Running 'scripts_seqModules/scripts_IonExpressSeq/project.single_IonExpressSeq.install_3.sh'" >> $logName;
echo "Variables passed via command-line from 'php/project.single_IonExpressSeq.install_2.php' :" >> $logName;
echo "\tuser     = '"$user"'" >> $logName;
echo "\tproject  = '"$project"'" >> $logName;
echo "\tmain_dir = '"$main_dir"'" >> $logName;
echo "#============================================================================== 3" >> $logName;

echo "#=====================================#" >> $logName;
echo "# Setting up locations and variables. #" >> $logName;
echo "#=====================================#" >> $logName;

echo "\tprojectDirectory = '$projectDirectory'" >> $logName;
echo "Setting up for processing." >> $condensedLog;

# Get genome and hapmap names in use from project's "genome.txt" file.
#    first line  => genome
#    second line => hapmap
genome=$(head -n 1 $projectDirectory"genome.txt");
hapmap=$(tail -n 1 $projectDirectory"genome.txt");
echo "\t'genome.txt' file entry." >> $logName;
echo "\t\tgenome = '"$genome"'" >> $logName;
if [ "$genome" = "$hapmap" ]
then
	hapmapInUse=0;
else
	echo "\t\thapmap = '"$hapmap"'" >> $logName;
	hapmapInUse=1;
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

# Get ploidy estimate from "ploidy.txt" in project directory.
ploidyEstimate=$(head -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyEstimate = '"$ploidyEstimate"'" >> $logName;

# Get ploidy baseline from "ploidy.txt" in project directory.
ploidyBase=$(tail -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyBase = '"$ploidyBase"'" >> $logName;

# Get parent name from "parent.txt" in project directory.
projectParent=$(head -n 1 $projectDirectory"parent.txt");
echo "\tparentProject = '"$projectParent"'" >> $logName;

# Define temporary directory for abra2 files.
abra2TempDirectory=$projectDirectory"abra2_temp/";
echo "\tabra2TempDirectory = '"$abra2TempDirectory"'" >> $logName;

echo "#============================================================================== 2" >> $logName;


if [ -f $projectDirectory"SNP_CNV_v1.txt" ]
then
	echo "\tDone: SAM -> BAM, new group headers, sorted." >> $logName;
	echo "\tBAM.indelrealignment done; Samtools.pileup generated." >> $logName;
else
	##==============================================================================
	## Trimming/cleanup of FASTQ files.
	##------------------------------------------------------------------------------
	echo "#=================================================#" >> $logName;
	echo "# Trimming of unbalanced FASTQ entries.           #" >> $logName;
	echo "#=================================================#" >> $logName;
	echo "Resolving FASTQ file errors." >> $condensedLog;
	currdir=$(pwd);
	cd $projectDirectory;
	sh $main_dir"scripts_seqModules/FASTQ_1_trimming.sh" $projectDirectory$datafile >> $logName;
	cd $currdir;
	echo "\tFASTQ files trimmed using : 'FASTQ_trimming.sh'" >> $logName;


	##==============================================================================
	## Initial processing of single-IonExpressSeq dataset.
	##------------------------------------------------------------------------------
	echo "#=========================================================#" >> $logName;
	echo "# Initial processing of single-end IonExpressSeq dataset. #" >> $logName;
	echo "#=========================================================#" >> $logName;

	# Align fastq against genome.
	echo "[[=- Align with Bowtie -=]]" >> $logName;
	echo "Aligning reads with Bowtie2 => SAM file." >> $condensedLog;

	threads=4;

	if [ -f $projectDirectory"data.bam" ]
	then
		echo "\tDone: SAM -> BAM, new group headers, sorted." >> $logName;
	else
		echo "\tBowtie : single-end reads aligning into SAM file." >> $logName;
		## Bowtie 2 command for single reads:
		echo "\nRunning bowtie2.\n";
		$bowtie2Directory"bowtie2" --very-sensitive -p $threads $genomeDirectory"bowtie_index" -U $projectDirectory$datafile -S $projectDirectory"data.sam";
			# -S : SAM output mode.
			# -p : number of threads to use.
			# -1 : dataset.
		    # --very-sensitive : a default set of configurations.
		echo "\tBowtie : single-end reads aligned into SAM file." >> $logName;

		echo "\tSamtools : converting Bowtie-SAM into compressed format (BAM) file." >> $logName;
		echo "Compressing SAM file => BAM file." >> $condensedLog;
		echo "\nRunning samtools:view.\n";
		$samtools_exec view -bT $genomeDirectory$genomeFASTA $projectDirectory"data.sam" > $projectDirectory"data.temp.bam";
		rm $projectDirectory"data.sam";
		echo "\tSamtools : Bowtie-SAM converted into compressed format (BAM) file." >> $logName;

		echo "\tPicard : Adding headers to Bowtie-BAM file." >> $logName;
		echo "Standardizing BAM read group headers." >> $condensedLog;
		currentDir=$(pwd);
		cd $picardDirectory;
		echo "\nRunning picard:AddOrReplaceReadGroups.\n";
		java -Xmx2g -jar AddOrReplaceReadGroups.jar INPUT=$projectDirectory"data.temp.bam" OUTPUT=$projectDirectory"data.bam" RGID=1 RGLB=1 RGPL=ILLUMINA RGPU=1 RGSM=SM;
		cd $currentDir;
		rm $projectDirectory"data.temp.bam";
		echo "\tPicard : Headers added to Bowtie-BAM file." >> $logName;

		echo "[[=- Sorting/Indexing BAM files -=]]" >> $logName;
		echo "\tSamtools : Bowtie-BAM sorting & indexing." >> $logName;
		echo "Sorting BAM file." >> $condensedLog;
		echo "\nRunning samtools:sort.\n";
		$samtools_exec sort $projectDirectory"data.bam" $projectDirectory"data_sorted";
		echo "Indexing BAM file." >> $condensedLog;
		echo "\nRunning samtools:index.\n";
		$samtools_exec index $projectDirectory"data_sorted.bam";
		echo "\tSamtools : Bowtie-BAM sorted & indexed." >> $logName;
	fi

	if [ -f $projectDirectory"data.pileup" ]
	then
		echo "\tBAM.indelrealignment done; Samtools.pileup generated." >> $logName;
	else
		#================================
		# Abra2: indel realignment.
		#--------------------------------
		echo "[[=- Indel realignment with ABRA2 analysis -=]]" >> $logName;
		echo "\tAbra2 : indel-realignment in process." >> $logName;
		echo "Indel realignment with ABRA2." >> $condensedLog;
		echo "\nRunning abra2.\n";
		ABRA2inputFile=$projectDirectory"data_sorted.bam";
		ABRA2outputFile=$projectDirectory"data_indelRealigned.bam";
		referenceFile=$genomeDirectory$genomeFASTA;
		mkdir $abra2TempDirectory;
		$java7Directory"java" -Xmx2g -jar $abra2_exec --in $ABRA2inputFile --out $ABRA2outputFile --ref $referenceFile --threads $cores --tmpdir $abra2TempDirectory > $projectDirectory"abra2.log";
		echo "\tAbra2 : indel-realignment done." >> $logName;
		# abra2-2.24.jar is missing file libAbra.so, which can be found in abra2-2.23.jar from github.com mozack/abra2.
		# example command-line from abra2 readme.
		# java -Xmx16G -jar abra2.jar --in input.bam --out output-sorted-realigned.bam --ref hg38.fa --threads 8 --targets targets.bed --tmpdir /your/tmpdir > abra.log
		# From paper: "Either the entire genome is traversed, or regions of interest can be specified via a bed file." in section 2.2.1 on page 2967.

		#================================
		# Sorting BAM file after Abra2.
		#--------------------------------
		echo "[[=- Sorting/Indexing BAM files -=]]" >> $logName;
		echo "\tSamtools : Bowtie-BAM sorting & indexing." >> $logName;
		echo "Sorting BAM file." >> $condensedLog;
		echo "\nRunning samtools:sort.\n";
		$samtools_exec sort -@ $cores $projectDirectory"data_indelRealigned.bam" -o $projectDirectory"data_sorted.bam" -T $projectDirectory;
		echo "Indexing BAM file." >> $condensedLog;
		echo "\nRunning samtools:index.\n";
		$samtools_exec index $projectDirectory"data_sorted.bam";
		echo "\tSamtools : Bowtie-BAM sorted & indexed." >> $logName;

		echo "#============================================================================== 3" >> $logName;

		echo "[[=- In-house SNP/CNV/INDEL analysis -=]]" >> $logName;
		usedFile=$projectDirectory"data_sorted.bam";
		echo "\tSamtools : Generating pileup.   (for SNP/CNV/INDEL analysis)" >> $logName;
		echo "Generating pileup file." >> $condensedLog;
		echo "\nRunning samtools:mpileup.\n";
		$samtools_exec mpileup -f $genomeDirectory$genomeFASTA $usedFile | awk '{print $1 " " $2 " " $3 " " $4 " " $5}' > $projectDirectory"data.pileup";
		echo "\tSamtools : Pileup generated." >> $logName;
	fi;

	echo "Processing pileup for CNVs, SNPs, & INDELs." >> $condensedLog;

	# ( echo "\tPython : Processing pileup for CNVs." >> $logName;
	# $python_exec $main_dir"py/counts_CNVs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_CNVs_v1.txt";
	# echo "\tPython : Pileup processed for CNVs." >> $logName; ) &
	#
	# ( echo "\tPython : Processing pileup for INDELs." >> $logName;
	# $python_exec $main_dir"py/counts_INDELs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_INDELS_v1.txt";
	# echo "\tPython : Pileup processed for INDELs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNPs." >> $logName;
	$python_exec $main_dir"scripts_seqModules/counts_SNPs_v5.py" $projectDirectory"data.pileup" > $projectDirectory"putative_SNPs_v4.txt";
	echo "\tPython : Pileup processed for SNPs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNP-CNV." >> $logName;
	$python_exec $main_dir"scripts_seqModules/counts_CNVs-SNPs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"SNP_CNV_v1.txt";
	echo "\tPython : Pileup processed for SNP-CNV." >> $logName; ) &

	wait;
fi

echo "Pileup processing is complete." >> $condensedLog;
echo "\nPileup processing complete.\n" >> $logName;
echo   "=========================================================================\n" >> $logName;

if [ $hapmapInUse = 0 ]
then
	echo "\nPassing processing on to 'project.IonExpressSeq.install_4.sh' for final analysis.\n" >> $logName;
	echo   "=================================================================================\n" >> $logName;
	sh $main_dir"scripts_seqModules/scripts_IonExpressSeq/project.IonExpressSeq.install_4.sh" $user $project $main_dir;
else
	echo "\nPassing processing on to 'project.IonExpressSeq.hapmap.install_4.sh' for final analysis.\n" >> $logName;
	echo   "========================================================================================\n" >> $logName;
	sh $main_dir"scripts_seqModules/scripts_IonExpressSeq/project.IonExpressSeq.hapmap.install_4.sh" $user $project $hapmap $main_dir;
fi
