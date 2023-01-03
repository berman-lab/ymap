#!/bin/bash -e
#
# project.single_WGseq.install_3.sh
#
set -e;
## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
project=$2;
main_dir=$(pwd)"/../../";


##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------

# Setup process_log.txt file.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";
echo "#.............................................................................." >> $logName;
echo "" >> $logName;
echo "Input to : project.single_WGseq.install_3.sh" >> $logName;
echo "\tuser     = "$user >> $logName;
echo "\tproject  = "$project >> $logName;
echo "\tmain_dir = "$main_dir >> $logName;
echo "" >> $logName;

# import locations of auxillary software for pipeline analysis.
. $main_dir"local_installed_programs.sh";
. $main_dir"config.sh";

# Define project directory.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";

echo "Running 'scripts_seqModules/scripts_WGseq/project.single_WGseq.install_3.sh'" >> $logName;
echo "Variables passed via command-line from 'scripts_seqModules/scripts_WGseq/project.single_WGseq.install_2.php' :" >> $logName;
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
echo "Location variables from 'genome.txt' file entry." >> $logName;
echo "\tgenome   = '"$genome"'" >> $logName;
if [ "$genome" = "$hapmap" ]
then
	hapmapInUse=0;
else
	echo "\thapmap   = '"$hapmap"'" >> $logName;
	hapmapInUse=1;
fi


if [ $hapmapInUse = 1 ]
then
	# Determine location of hapmap being used.
	if [ -d $main_dir"users/"$user"/hapmaps/"$hapmap"/" ]
	then
		hapmapDirectory=$main_dir"users/"$user"/hapmaps/"$hapmap"/";
		hapmapUser=$user;
	elif [ -d $main_dir"users/default/hapmaps/"$hapmap"/" ]
	then
		hapmapDirectory=$main_dir"users/default/hapmaps/"$hapmap"/";
		hapmapUser="default";
	fi
	echo "\thapmapDirectory = '"$hapmapDirectory"'" >> $logName;
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

## Define temporary directory for FASTQC files.
#fastqcTempDirectory=$projectDirectory"fastqc_temp/";
#echo "\tfastqcTempDirectory = '"$fastqcTempDirectory"'" >> $logName;

# Get ploidy estimate from "ploidy.txt" in project directory.
ploidyEstimate=$(head -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyEstimate = '"$ploidyEstimate"'" >> $logName;

# Get ploidy baseline from "ploidy.txt" in project directory.
ploidyBase=$(tail -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyBase = '"$ploidyBase"'" >> $logName;

# Get parent name from "parent.txt" in project directory.
projectParent=$(head -n 1 $projectDirectory"parent.txt");
echo "\tparentProject = '"$projectParent"'" >> $logName;

# Define temporary directory for abra2C files.
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
	echo "#=======================================================================================#" >> $logName;
	echo "# Trimming of unbalanced FASTQ entries using 'scripts_seqModules/FASTQ_1_trimming.sh'.  #" >> $logName;
	echo "#=======================================================================================#" >> $logName;
	echo "Resolving FASTQ file errors." >> $condensedLog;
	currdir=$(pwd);
	cd $projectDirectory;
	sh $main_dir"scripts_seqModules/FASTQ_1_trimming.sh" $projectDirectory$datafile >> $logName;
	cd $currdir;

	##==============================================================================
	## Initial processing of single-WGseq dataset.
	##------------------------------------------------------------------------------
	echo "#=================================================#" >> $logName;
	echo "# Initial processing of single-end WGseq dataset. #" >> $logName;
	echo "#=================================================#" >> $logName;

	# Align fastq against genome.
	echo "[[=- Align with Bowtie -=]]" >> $logName;
	echo "Aligning reads with Bowtie2 => SAM file." >> $condensedLog;

	if [ -f $projectDirectory"data.bam" ]
	then
		echo "\tDone: SAM -> BAM, new group headers, sorted." >> $logName;
	else
		echo "\tBowtie : single-end reads aligning into SAM file." >> $logName;
		## Bowtie 2 command for single reads:
		echo "\nRunning bowtie2.\n";
		echo "\tbowtie2 --very-sensitive -p "$cores" "$genomeDirectory"bowtie_index -U "$projectDirectory$datafile" -S "$projectDirectory"data.sam;" >> $logName;
		$bowtie2Directory"bowtie2" --very-sensitive -p $cores $genomeDirectory"bowtie_index" -U $projectDirectory$datafile -S $projectDirectory"data.sam";
			# -S : SAM output mode.
			# -p : number of threads to use.
			# -1 : dataset.
		    # --very-sensitive : a default set of configurations.
		echo "\tBowtie : single-end reads aligned into SAM file." >> $logName;

		echo "\tSamtools : converting Bowtie-SAM into compressed format (BAM) file." >> $logName;
		echo "Compressing SAM file => BAM file." >> $condensedLog;
		echo "\nRunning samtools:view.\n";
		$samtools_exec view -@ $cores -bT $genomeDirectory$genomeFASTA $projectDirectory"data.sam" > $projectDirectory"data.temp.bam";
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
		$samtools_exec sort -@ $cores $projectDirectory"data.bam" -o $projectDirectory"data_sorted.bam" -T $projectDirectory;
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
# Removed FastQC and GATK.
#--------------------------------
#		echo "[[=- GATK analysis, indel-realignment -=]]" >> $logName;
#		GATKinputFile=$projectDirectory"data_sorted.bam";
#		GATKoutputFile1=$projectDirectory"data_forIndelRealigner.intervals";
#		GATKoutputFile2=$projectDirectory"data_indelRealigned.bam";
#		GATKreference=$genomeDirectory$genomeFASTA;
#
#		# 'LENIENT' should allow GATK to ignore problem reads.
#		GATKoptions="-S LENIENT -filterMBQ";
#
#		## FASTQC : determine read quality format in use.
#		echo "\tFASTQC : Read quality coding is required." >> $logName;
#		## FASTQC analysis of FASTQ input files to determine read quality coding.
#		## This is needed for determining GATK command option to deal with alternate coding.
#		#		S - Sanger        Phred+33,  raw reads typically (0, 40)
#		#		X - Solexa        Solexa+64, raw reads typically (-5, 40)
#		#		I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
#		#		J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
#		#		L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
#		echo "\tWorking directory : '"$fastqcTempDirectory"'" >> $logName;
#
#		mkdir $fastqcTempDirectory;
#
#		# -o : points to standardized ouptut temp directory.
#		echo "Identifying quality coding used in FASTA." >> $condensedLog;
#		FASTQClog1=$projectDirectory"fastqc.process.log";
#		echo "\nRunning fastqc.\n";
#               echo "\tFASTQC : "$fastqcDirectory"fastqc -t "$cores" -o "$fastqcTempDirectory" "$projectDirectory$datafile" > "$FASTQClog1";" >> $logName;
#		$fastqcDirectory"fastqc" -t $cores -o $fastqcTempDirectory $projectDirectory$datafile > $FASTQClog1;
#		echo "\tFASTQC log output :" >> $logName;
#		sed 's/^/\t\t|/;' $FASTQClog1 >> $logName;
#		rm $FASTQClog1;
#
#		# get file extension.
#		fileExt=${datafile#*.};
#
#		# generate file name with extension removed and '_fastqc' appended, to match FASTQC output files.
#		findStr="."$fileExt;
#		replaceStr="_fastqc";
#		FASTQC_temp_file=$(echo $datafile | sed -e "s/$findStr/$replaceStr/g");
#
#		# parse read quality encoding format from FASTQC output.
#		qualityCodingLine=`sed -n 6,6'p' $fastqcTempDirectory$FASTQC_temp_file"/fastqc_data.txt"`;
#		qualityCoding=$(echo $qualityCodingLine | cut -c9-);
#
#		# cleanup extraneous FASTQC files.
#		rm $fastqcTempDirectory$FASTQC_temp_file".zip";
#		rm -rf $fastqcTempDirectory$FASTQC_temp_file;
#
#		# Trim starting and ending whitespace from quality coding string.
#		trimmedQualityCoding="${qualityCoding#"${qualityCoding%%[![:space:]]*}"}";                # remove leading whitespace characters
#		trimmedQualityCoding="${trimmedQualityCoding%"${trimmedQualityCoding##*[![:space:]]}"}";  # remove trailing whitespace characters
#		echo "\t\tQuality coding method = '"$trimmedQualityCoding"'." >> $logName;
#
#		# Expected output:
#		#       "Illumina 1.5"          : starting with 64.
#		#       "Sanger / Illumina 1.9" : normal, starting with 33.   [MiSeq]
#
#		# If quality coding string matches a type known to start with 64, update GATK option string.
#		if [ "$trimmedQualityCoding" = "Illumina 1.5" ]
#		then
#			GATKoptions=$GATKoptions" -fixMisencodedQuals";
#			echo "\t\tGATK option to deal with quality coding = ' -fixMisencodedQuals'." >> $logName;
#		fi
#
#		GATKlog1=$projectDirectory"gatk.RealignerTargetCreator.log";
#		GATKlog2=$projectDirectory"gatk.IndelRealigner.log";
#
#		echo "\tGATK options = '"$GATKoptions"'" >> $logName;
#		echo "\tGATK : preparing for IndelRealignment." >> $logName;
#		echo "Preparing for indel realignment." >> $condensedLog;
#		echo "\nRunning gatk:RealignerTargetCreator.\n";
#		$java7Directory"java" -Xmx2g -jar $gatkDirectory"GenomeAnalysisTK.jar" -T RealignerTargetCreator -nt $cores -I $GATKinputFile -R $GATKreference -o $GATKoutputFile1 $GATKoptions > $GATKlog1;
#		sed 's/^/\t\t|/;' $GATKlog1 >> $logName;
#		echo "\tGATK : prepared for IndelRealignment." >> $logName;
#		echo "\tGATK : performing IndelRealignment." >> $logName;
#		echo "Realigning indels." >> $condensedLog;
#		echo "\nRunning gatk:IndelRealigner.\n";
#		$java7Directory"java" -Xmx2g -jar $gatkDirectory"GenomeAnalysisTK.jar" -T IndelRealigner -I $GATKinputFile -R $GATKreference -targetIntervals $GATKoutputFile1 -o $GATKoutputFile2 $GATKoptions > $GATKlog2;
#		sed 's/^/\t\t|/;' $GATKlog2 >> $logName;
#		echo "\tGATK : performed IndelRealignment." >> $logName;
#
#		rm $GATKlog1;
#		rm $GATKlog2;

                #================================
                # Abra2: indel realignment.
                #--------------------------------
                echo "[[=- ABRA2 analysis, indel-realignment -=]]" >> $logName;
                echo "Indel realignment with ABRA2." >> $condensedLog;
                echo "\nRunning abra2.\n";
                ABRA2inputFile=$projectDirectory"data_sorted.bam";
                ABRA2outputFile=$projectDirectory"data_indelRealigned.bam";
                referenceFile=$genomeDirectory$genomeFASTA;
                mkdir $abra2TempDirectory;
                echo "__";
                echo $java7Directory"java -Xmx2g -jar "$abra2_exec" --in "$ABRA2inputFile" --out "$ABRA2outputFile" --ref "$referenceFile" --threads "$cores" --tmpdir "$abra2TempDirectory" > "$projectDirectory"abra2.log" >> $logName;
                echo "__";
                $java7Directory"java" -Xmx2g -jar $abra2_exec --in $ABRA2inputFile --out $ABRA2outputFile --ref $referenceFile --threads $cores --tmpdir $abra2TempDirectory > $projectDirectory"abra2.log";
                # abra2-2.24.jar is missing file libAbra.so, which can be found in abra2-2.23.jar from github.com mozack/abra2.
                # example command-line from abra2 readme.
                # java -Xmx16G -jar abra2.jar --in input.bam --out output-sorted-realigned.bam --ref hg38.fa --threads 8 --targets targets.bed --tmpdir /your/tmpdir > abra.log
                # From paper: "Either the entire genome is traversed, or regions of interest can be specified via a bed file." in section 2.2.1 on page 2967.

		echo "#============================================================================== 3" >> $logName;

		echo "[[=- In-house SNP/CNV/INDEL analysis -=]]" >> $logName;

                # data_sorted.bam is to be used if no indel-realignment is done. data_indelRealigned.bam if indel-realignment done.
		usedFile=$projectDirectory"data_indelRealigned.bam";
#               usedFile=$projectDirectory"data_sorted.bam";

		echo "\tSamtools : Generating pileup.   (for SNP/CNV/INDEL analysis)" >> $logName;
		echo "Generating pileup file." >> $condensedLog;
		echo "\nRunning samtools:mpileup.\n";
		$python_exec $main_dir"scripts_seqModules/parallel_mpileup.py" $samtools_exec $genomeDirectory$genomeFASTA $usedFile $logName $cores $genomeDirectory data.pileup 2>> $logName;
		echo "\tSamtools : Pileup generated." >> $logName;
		echo "Processing pileup for CNVs, SNPs, & INDELs." >> $condensedLog;
	fi

	# ( echo "\tPython : Processing pileup for CNVs." >> $logName;
	# $python_exec $main_dir"scripts_seqModules/counts_CNVs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_CNVs_v1.txt";
	# echo "\tPython : Pileup processed for CNVs." >> $logName; ) &
	#
	# ( echo "\tPython : Processing pileup for INDELs." >> $logName;
	# $python_exec $main_dir"scripts_seqModules/counts_INDELs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_INDELS_v1.txt";
	# echo "\tPython : Pileup processed for INDELs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNPs." >> $logName;
	$python_exec $main_dir"scripts_seqModules/counts_SNPs_v5.py" $projectDirectory"data.pileup" > $projectDirectory"putative_SNPs_v4.txt";
	echo "\tPython : Pileup processed for SNPs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNP-CNV." >> $logName;
	$python_exec $main_dir"scripts_seqModules/counts_CNVs-SNPs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"SNP_CNV_v1.txt";
	echo "\tPython : Pileup processed for SNP-CNV." >> $logName; ) &

	wait;
fi
if [ $hapmapInUse = 1 ]
then
	if [ -f $projectDirectory"trimmed_SNPs_v5.txt" ]
	then
		echo "\tPython : Simplify child putative_SNP list to contain only those loci found in the haplotype map." >> $logName;
		echo "\t\tDone." >> $logName;
	else
		echo "\tPython : Simplify child putative_SNP list to contain only those loci found in the haplotype map." >> $logName;
		echo "\t\t| Inputs to python script:" >> $logName;
		echo "\t\t|\tgenome     = '$genome'"     >> $logName;
		echo "\t\t|\tgenomeUser = '$genomeUser'" >> $logName;
		echo "\t\t|\tproject    = '$project'"    >> $logName;
		echo "\t\t|\tuser       = '$user'"       >> $logName;
		echo "\t\t|\thapmap     = '$hapmap'"     >> $logName;
		echo "\t\t|\thapmapUser = '$hapmapUser'" >> $logName;
		echo "\t\t|\tmain_dir   = '$main_dir'"   >> $logName;
		$python_exec $main_dir"scripts_seqModules/putative_SNPs_from_hapmap_in_child.py" $genome $genomeUser $project $user $hapmap $hapmapUser $main_dir > $projectDirectory"trimmed_SNPs_v5.txt"
		echo "\t\tDone." >> $logName;
	fi
fi

echo "Pileup processing is complete." >> $condensedLog;
echo "\nPileup processing complete.\n" >> $logName;
echo   "=========================================================================\n" >> $logName;

if [ $hapmapInUse = 0 ]
then
	echo "Passing processing on to 'scripts_seqModules/scripts_WGseq/project.WGseq.install_4.sh' for final analysis." >> $logName;
	echo "\t"$main_dir"scripts_seqModules/scripts_WGseq/project.WGseq.install_4.sh "$user" "$project" "$main_dir >> $logName;
	echo "Script executed from: ${PWD}" >> $logName;

	echo   "=========================================================================\n" >> $logName;
	sh $main_dir"scripts_seqModules/scripts_WGseq/project.WGseq.install_4.sh" $user $project $main_dir;
else
	echo "Passing processing on to 'scripts_seqModules/scripts_WGseq/project.WGseq.hapmap.install_4.sh' for final analysis." >> $logName;
	echo "\t"$main_dir"scripts_seqModules/scripts_WGseq/project.WGseq.hapmap.install_4.sh "$user" "$project" "$hapmap" "$main_dir >> $logName;
	echo "Script executed from: ${PWD}" >> $logName;

	echo "=========================================================================\n" >> $logName;
	sh $main_dir"scripts_seqModules/scripts_WGseq/project.WGseq.hapmap.install_4.sh" $user $project $hapmap $main_dir;
fi
