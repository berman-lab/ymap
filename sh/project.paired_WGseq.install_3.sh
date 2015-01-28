#!/bin/bash -e
#
# project.paired_WGseq.install_3.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
main_dir=$(pwd)"/../";

#user="darren";
#project="test_12353_vs_hapmap_2b";
#main_dir="/heap/hapmap/bermanlab/";

#user="morlurie";
#project="26_try4";
#main_dir="/heap/hapmap/bermanlab/";

##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------

# import locations of auxillary software for pipeline analysis.
. $main_dir"sh/local_installed_programs.sh";

# Define project directory.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";

# Setup process_log.txt file.
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";
#chmod 0755 $logName;
echo "#.............................................................................." >> $logName;
echo "Running 'sh/project.paired_WGseq.install_3.sh'" >> $logName;
echo "Variables passed via command-line from 'php/project.paired_WGseq.install_2.php' :" >> $logName;
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

# Get first data file name from "datafiles.txt";
datafile1=$(head -n 1 $projectDirectory"datafiles.txt");
echo "\tdatafile 1 = '"$datafile1"'" >> $logName;

# Get second data file name from "datafiles.txt";
datafile2=$(tail -n 1 $projectDirectory"datafiles.txt");
echo "\tdatafile 2 = '"$datafile2"'" >> $logName;

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
	sh $main_dir"sh/FASTQ_2_trimming.sh" $projectDirectory$datafile1 $projectDirectory$datafile2 >> $logName;
	cd $currdir;
	echo "\tFASTQ files trimmed using : 'FASTQ_trimming.sh'" >> $logName;


	##==============================================================================
	## Initial processing of paired-WGseq dataset.
	##------------------------------------------------------------------------------
	echo "#=================================================#" >> $logName;
	echo "# Initial processing of paired-end WGseq dataset. #" >> $logName;
	echo "#=================================================#" >> $logName;

	# Align fastq against genome.
	echo "[[=- Align with Bowtie -=]]" >> $logName;
	echo "Aligning reads with Bowtie2 => SAM file." >> $condensedLog;

	threads=4;

	if [ -f $projectDirectory"data.bam" ]
	then
		echo "\tDone: SAM -> BAM, new group headers, sorted." >> $logName;
	else
		echo "\tBowtie : paired-end reads aligning into SAM file." >> $logName;
		## Bowtie 2 command for paired reads:
		echo "\nRunning bowtie2.\n";
		bowtie2 --very-sensitive -p $threads $genomeDirectory"bowtie_index" -1 $projectDirectory$datafile1 -2 $projectDirectory$datafile2 -S $projectDirectory"data.sam";
			# -S : SAM output mode.
			# -p : number of threads to use.
			# -1 : dataset.
		    # --very-sensitive : a default set of configurations.
		echo "\tBowtie : paired-end reads aligned into SAM file." >> $logName;

		echo "\tSamtools : converting Bowtie-SAM into compressed format (BAM) file." >> $logName;
		echo "Compressing SAM file => BAM file." >> $condensedLog;
		echo "\nRunning samtools:view.\n";
		samtools view -bT $genomeDirectory$genomeFASTA $projectDirectory"data.sam" > $projectDirectory"data.temp.bam";
		rm $projectDirectory"data.sam";
		echo "\tSamtools : Bowtie-SAM converted into compressed format (BAM) file." >> $logName;

		echo "\tPicard : Adding headers to Bowtie-BAM file." >> $logName;
		echo "Standardizing BAM read group headers." >> $condensedLog;
		echo "\nRunning picard:AddOrReplaceReadGroups.\n";
		java -Xmx2g -jar $picardDirectory"AddOrReplaceReadGroups.jar" INPUT=$projectDirectory"data.temp.bam" OUTPUT=$projectDirectory"data.bam" RGID=1 RGLB=1 RGPL=ILLUMINA RGPU=1 RGSM=SM;
		rm $projectDirectory"data.temp.bam";
		echo "\tPicard : Headers added to Bowtie-BAM file." >> $logName;

		echo "[[=- Sorting/Indexing BAM files -=]]" >> $logName;
		echo "\tSamtools : Bowtie-BAM sorting & indexing." >> $logName;
		echo "Sorting BAM file." >> $condensedLog;
		echo "\nRunning samtools:sort.\n";
		samtools sort $projectDirectory"data.bam" $projectDirectory"data_sorted";
		echo "Indexing BAM file." >> $condensedLog;
		echo "\nRunning samtools:index.\n";
		samtools index $projectDirectory"data_sorted.bam";
		echo "\tSamtools : Bowtie-BAM sorted & indexed." >> $logName;
	fi;

	if [ -f $projectDirectory"data.pileup" ]
	then
		echo "\tBAM.indelrealignment done; Samtools.pileup generated.." >> $logName;
	else
		echo "[[=- GATK analysis, indel-realignment -=]]" >> $logName;
		GATKinputFile=$projectDirectory"data_sorted.bam";
		GATKoutputFile1=$projectDirectory"data_forIndelRealigner.intervals";
		GATKoutputFile2=$projectDirectory"data_indelRealigned.bam";
		GATKreference=$genomeDirectory$genomeFASTA;

		# 'LENIENT' should allow GATK to ignore problem reads.
		GATKoptions="-S LENIENT -filterMBQ";

		## FASTQC : determine read quality format in use.
		echo "\tFASTQC : Read quality coding is required." >> $logName;
		## FASTQC analysis of FASTQ input files to determine read quality coding.
		## This is needed for determining GATK command option to deal with alternate coding.
		#		S - Sanger        Phred+33,  raw reads typically (0, 40)
		#		X - Solexa        Solexa+64, raw reads typically (-5, 40)
		#		I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
		#		J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
		#		L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
		echo "\tWorking directory : '"$fastqcTempDirectory"'" >> $logName;

		mkdir $fastqcTempDirectory;

		# -o : points to standardized ouptut temp directory.
		echo "Identifying quality coding used in FASTA." >> $condensedLog;
		FASTQClog1=$projectDirectory"fastqc.process.log";
		echo "\nRunning fastqc.\n";
		$fastqcDirectory"fastqc" -o $fastqcTempDirectory $projectDirectory$datafile1 > $FASTQClog1;
		echo "\tFASTQC log output :" >> $logName;
		sed 's/^/\t\t|/;' $FASTQClog1 >> $logName;
		rm $FASTQClog1;

		echo "fastQC_1" >> $logName;

		# get file extension.
		fileExt=${datafile1#*.};

		echo "fastQC_2" >> $logName;

		# generate file name with extension removed and '_fastqc' appended, to match FASTQC output files.
		findStr="."$fileExt;
		replaceStr="_fastqc";
		FASTQC_temp_file=$(echo $datafile1 | sed -e "s/$findStr/$replaceStr/g");

		echo "fastQC_3" >> $logName;

		# parse read quality encoding format from FASTQC output.
		qualityCodingLine=`sed -n 6,6'p' $fastqcTempDirectory$FASTQC_temp_file"/fastqc_data.txt"`;
		qualityCoding=$(echo $qualityCodingLine | cut -c9-);

		echo "fastQC_4" >> $logName;

		# cleanup extraneous FASTQC files.
		rm $fastqcTempDirectory$FASTQC_temp_file".zip";
		rm -rf $fastqcTempDirectory$FASTQC_temp_file;

		echo "fastQC_5" >> $logName;

		# Trim starting and ending whitespace from quality coding string.
		trimmedQualityCoding="${qualityCoding#"${qualityCoding%%[![:space:]]*}"}";                # remove leading whitespace characters
		trimmedQualityCoding="${trimmedQualityCoding%"${trimmedQualityCoding##*[![:space:]]}"}";  # remove trailing whitespace characters
		echo "\t\tQuality coding method = '"$trimmedQualityCoding"'." >> $logName;

		# Expected output:
		#       "Illumina 1.5"          : starting with 64.
		#       "Sanger / Illumina 1.9" : normal, starting with 33.   [MiSeq]

		# If quality coding string matches a type known to start with 64, update GATK option string.
		if [ "$trimmedQualityCoding" = "Illumina 1.5" ]
		then
			GATKoptions=$GATKoptions" -fixMisencodedQuals";
			echo "\t\tGATK option to deal with quality coding = ' -fixMisencodedQuals'." >> $logName;
		fi

		GATKlog1=$projectDirectory"gatk.RealignerTargetCreator.log";
		GATKlog2=$projectDirectory"gatk.IndelRealigner.log";
		currentDir=$(pwd);
		cd $gatkDirectory;

		echo "\tGATK options = '"$GATKoptions"'" >> $logName;
		echo "\tGATK : preparing for IndelRealignment." >> $logName;
		echo "Preparing for indel realignment." >> $condensedLog;
		echo "\nRunning gatk:RealignerTargetCreator.\n";
		$java7Directory"java" -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -I $GATKinputFile -R $GATKreference -o $GATKoutputFile1 $GATKoptions > $GATKlog1;
		sed 's/^/\t\t|/;' $GATKlog1 >> $logName;
		echo "\tGATK : prepared for IndelRealignment." >> $logName;
		echo "\tGATK : performing IndelRealignment." >> $logName;
		echo "Realigning indels." >> $condensedLog;
		echo "\nRunning gatk:IndelRealigner.\n";
		$java7Directory"java" -jar GenomeAnalysisTK.jar -T IndelRealigner -I $GATKinputFile -R $GATKreference -targetIntervals $GATKoutputFile1 -o $GATKoutputFile2 $GATKoptions > $GATKlog2;
		sed 's/^/\t\t|/;' $GATKlog2 >> $logName;
		echo "\tGATK : performed IndelRealignment." >> $logName;
		rm $GATKlog1;
		rm $GATKlog2;

		cd $currentDir;

		echo "#============================================================================== 3" >> $logName;

		echo "[[=- In-house SNP/CNV/INDEL analysis -=]]" >> $logName;
		usedFile=$projectDirectory"data_indelRealigned.bam";

		echo "\tSamtools : Generating pileup.   (for SNP/CNV/INDEL analysis)" >> $logName;
		echo "Generating pileup file." >> $condensedLog;
		echo "\nRunning samtools:mpileup.\n";
		samtools mpileup -f $genomeDirectory$genomeFASTA $usedFile | awk '{print $1 " " $2 " " $3 " " $4 " " $5}' > $projectDirectory"data.pileup";
		echo "\tSamtools : Pileup generated." >> $logName;
	fi;

	echo "Processing pileup for CNVs & SNPs." >> $condensedLog;

	# ( echo "\tPython : Processing pileup for CNVs." >> $logName;
	# $python_exec $main_dir"py/counts_CNVs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_CNVs_v1.txt";
	# echo "\tPython : Pileup processed for CNVs." >> $logName; ) &
	#
	# ( echo "\tPython : Processing pileup for INDELs." >> $logName;
	# $python_exec $main_dir"py/counts_INDELs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_INDELS_v1.txt";
	# echo "\tPython : Pileup processed for INDELs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNPs." >> $logName;
	$python_exec $main_dir"py/counts_SNPs_v5.py" $projectDirectory"data.pileup" > $projectDirectory"putative_SNPs_v4.txt";
	echo "\tPython : Pileup processed for SNPs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNP-CNV." >> $logName;
	$python_exec $main_dir"py/counts_CNVs-SNPs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"SNP_CNV_v1.txt";
	echo "\tPython : Pileup processed for SNP-CNV." >> $logName; ) &

	wait;
fi

echo "Pileup processing is complete." >> $condensedLog;
echo "\nPileup processing complete.\n" >> $logName;

if [ $hapmapInUse = 0 ]
then
	echo "\nPassing processing on to 'project.WGseq.install_4.sh' for final analysis.\n" >> $logName;
	echo   "=========================================================================\n" >> $logName;
	sh $main_dir"sh/project.WGseq.install_4.sh" $user $project $main_dir;
else
	echo "\nPassing processing on to 'project.WGseq.hapmap.install_4.sh' for final analysis.\n" >> $logName;
	echo   "================================================================================\n" >> $logName;
	sh $main_dir"sh/project.WGseq.hapmap.install_4.sh" $user $project $hapmap $main_dir;
fi
