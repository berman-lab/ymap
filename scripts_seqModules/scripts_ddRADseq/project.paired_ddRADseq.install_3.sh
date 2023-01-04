#!/bin/bash -e
#
# project.paired_ddRADseq.install_3.sh
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

# import locations of auxillary software for pipeline analysis.
. $main_dir"local_installed_programs.sh";
. $main_dir"config.sh";

# Define project directory.
projectDirectory=$main_dir"users/"$user"/projects/"$project"/";

# Setup process_log.txt file.
logName=$projectDirectory"process_log.txt";
condensedLog=$projectDirectory"condensed_log.txt";
chmod 0755 $logName;
echo "#.............................................................................." >> $logName;
echo "Running 'scripts_seqModules/scripts_ddRADseq/project.paired_ddRADseq.install_3.sh'" >> $logName;
echo "Variables passed via command-line from 'scripts_seqModules/scripts_ddRADseq/project.paired_ddRADseq.install_2.php' :" >> $logName;
echo "\tuser     = '"$user"'" >> $logName;
echo "\tproject  = '"$project"'" >> $logName;
echo "\tmain_dir = '"$main_dir"'" >> $logName;
echo "#============================================================================== 3" >> $logName;

echo "#=====================================#" >> $logName;
echo "# Setting up locations and variables. #" >> $logName;
echo "#=====================================#" >> $logName;

echo "\tprojectDirectory = '$projectDirectory'" >> $logName;
echo "Setting up for processing." >> $condensedLog;

# Get setup information from project files.
# "genome.txt"
#    first line  => genome
#    second line => hapmap
# "dataType.txt"
#    5th character, 0=no indel-realignment, 1= indel-realignment.
genome=$(head -n 1 $projectDirectory"genome.txt");
hapmap=$(tail -n 1 $projectDirectory"genome.txt");
dataType=$(head -n 1 $projectDirectory"dataType.txt");
echo "\t'genome.txt' file entry." >> $logName;
echo "\t\tgenome = '"$genome"'" >> $logName;
if [ "$genome" = "$hapmap" ]
then
	hapmapInUse=0;
else
	echo "\t\thapmap = '"$hapmap"'" >> $logName;
	hapmapInUse=1;
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
fi
indelrealign_bool=$(echo $dataType | cut -c5-5);  # 0=no indel-realignment; 1=indel-realignment.

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

# Get ploidy estimate from "ploidy.txt" in project directory.
ploidyEstimate=$(head -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyEstimate = '"$ploidyEstimate"'" >> $logName;

# Get ploidy baseline from "ploidy.txt" in project directory.
ploidyBase=$(tail -n 1 $projectDirectory"ploidy.txt");
echo "\tploidyBase = '"$ploidyBase"'" >> $logName;

# Get parent name from "parent.txt" in project directory.
projectParent=$(head -n 1 $projectDirectory"parent.txt");
echo "\tparentProject = '"$projectParent"'" >> $logName;

if [ $indelrealign_bool = 1 ]
then
	# Define temporary directory for abra2 files.
	abra2TempDirectory=$projectDirectory"abra2_temp/";
	echo "\tabra2TempDirectory = '"$abra2TempDirectory"'" >> $logName;
fi

# Determine location of parent being used.
if [ -d $main_dir"users/"$user"/projects/"$projectParent"/" ]
then
	projectParentDirectory=$main_dir"users/"$user"/genomes/"$projectParent"/";
	projectParentUser=$user;
elif [ -d $main_dir"users/default/projects/"$projectParent"/" ]
then
	projectParentDirectory=$main_dir"users/default/genomes/"$projectParent"/";
	projectParentUser="default";
fi


reflocation=$main_dir"users/"$genomeUser"/genomes/"$genome"/";                 # Directory where FASTA file is kept.
FASTA=`sed -n 1,1'p' $reflocation"reference.txt"`;                             # Name of FASTA file.
FASTAname=$(echo $FASTA | sed 's/.fasta//g');                                  # name of genome file, without file type.
RestrctionEnzymes=`sed -n 1,1'p' $projectDirectory"restrictionEnzymes.txt"`;   # Name of restriction enxyme list file.
ddRADseq_FASTA=$FASTAname"."$RestrctionEnzymes".fasta";                        # Name of digested reference for ddRADseq analysis, using chosen restriction enzymes.


echo "#============================================================================== 2" >> $logName;


if [ -f $projectDirectory"SNP_CNV_v1.txt" ]
then
	echo "\tDone: SAM -> BAM, new group headers, sorted." >> $logName;
	echo "\tDone: BAM.indelrealignment." >> $logName;
	echo "\tDone: Samtools.pileup." >> $logName;
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
	sh $main_dir"scripts_seqModules/FASTQ_2_trimming.sh" $projectDirectory$datafile1 $projectDirectory$datafile2 >> $logName;
	cd $currdir;
	echo "\tFASTQ files trimmed using : 'FASTQ_trimming.sh'" >> $logName;


	##==============================================================================
	## Initial processing of paired-ddRADseq dataset.
	##------------------------------------------------------------------------------
	echo "#====================================================#" >> $logName;
	echo "# Initial processing of paired-end ddRADseq dataset. #" >> $logName;
	echo "#====================================================#" >> $logName;

	# Align fastq against genome.
	echo "[[=- Align with Bowtie -=]]" >> $logName;
	echo "Aligning reads with Bowtie2 => SAM file." >> $condensedLog;

	if [ -f $projectDirectory"data.bam" ]
	then
		echo "\tDone: SAM -> BAM, new group headers, sorted." >> $logName;
	else
		echo "\tBowtie : paired-end reads aligning into SAM file." >> $logName;
		## Bowtie 2 command for paired reads:
		echo "\nRunning bowtie2.\n";
		$bowtie2Directory"bowtie2" --very-sensitive -p $cores $genomeDirectory"bowtie_index" -1 $projectDirectory$datafile1 -2 $projectDirectory$datafile2 -S $projectDirectory"data.sam";
			# -S : SAM output mode.
			# -p : number of threads to use.
			# -1 : dataset.
		    # --very-sensitive : a default set of configurations.
		echo "\tBowtie : paired-end reads aligned into SAM file." >> $logName;

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
		echo "\tBAM.indelrealignment done; Samtools.pileup generated.." >> $logName;
	else
		if [ $indelrealign_bool = 1 ]
		then
			#================================
			# Abra2: indel realignment.
			#--------------------------------
			echo "[[=- Indel realignment with ABRA2 -=]]" >> $logName;
			echo "\tAbra2 : indel-realignment in process." >> $logName;
			echo "Indel realignment with ABRA2." >> $condensedLog;
			echo "\nRunning abra2.\n";
			ABRA2bedFile=$genomeDirectory"genome.bed";
			ABRA2inputFile=$projectDirectory"data_sorted.bam";
			ABRA2outputFile=$projectDirectory"data_indelRealigned.bam";
			referenceFile=$genomeDirectory$genomeFASTA;
			mkdir $abra2TempDirectory;
			echo ""  >> $logName;
			echo "command: "$java7Directory"java -Xmx16g -jar "$abra2_exec" --in "$ABRA2inputFile" --out "$ABRA2outputFile" --ref "$referenceFile" --threads "$cores" --targets "$ABRA2bedFile" --tmpdir "$abra2TempDirectory" > "$projectDirectory"abra2.log"  >> $logName;
			echo ""  >> $logName;
			$java7Directory"java" -Xmx16g -jar $abra2_exec --in $ABRA2inputFile --out $ABRA2outputFile --ref $referenceFile --threads $cores --targets $ABRA2bedFile --tmpdir $abra2TempDirectory > $projectDirectory"abra2.log";
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
		else
			echo "[[=- Indel realignment not being done -=]]" >> $logName;
		fi

		echo "#============================================================================== 3" >> $logName;

		echo "[[=- In-house SNP/CNV/INDEL analysis -=]]" >> $logName;
		usedFile=$projectDirectory"data_sorted.bam";
		echo "\tSamtools : Generating pileup.   (for SNP/CNV/INDEL analysis)" >> $logName;
		echo "Generating pileup file." >> $condensedLog;
		echo "\nRunning samtools:mpileup.\n";
		$python_exec $main_dir"scripts_seqModules/parallel_mpileup.py" $samtools_exec $genomeDirectory$genomeFASTA $usedFile $logName $cores $genomeDirectory data.pileup 2>> $logName;
		echo "\tSamtools : Pileup generated." >> $logName;
	fi

	echo "Processing pileup for CNVs & SNPs." >> $condensedLog;

	# ( echo "\tPython : Processing pileup for CNVs." >> $logName;
	# $python_exec $main_dir"scripts_seqModules/scripts_ddRADseq/counts_CNVs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_CNVs_v1.txt";
	# echo "\tPython : Pileup processed for CNVs." >> $logName; ) &
	#
	# ( echo "\tPython : Processing pileup for INDELs." >> $logName;
	# $python_exec $main_dir"scripts_seqModules/scripts_ddRADseq/counts_INDELs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"putative_INDELS_v1.txt";
	# echo "\tPython : Pileup processed for INDELs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNPs." >> $logName;
	$python_exec $main_dir"scripts_seqModules/counts_SNPs_v5.py" $projectDirectory"data.pileup" > $projectDirectory"putative_SNPs_v4.txt";
	echo "\tPython : Pileup processed for SNPs." >> $logName; ) &

	( echo "\tPython : Processing pileup for SNP-CNV." >> $logName;
	$python_exec $main_dir"scripts_seqModules/counts_CNVs-SNPs_v1.py" $projectDirectory"data.pileup" > $projectDirectory"SNP_CNV_v1.txt";
	echo "\tPython : Pileup processed for SNP-CNV." >> $logName; ) &

	wait;
fi
if [ -f $projectDirectory"trimmed_SNPs_v4.txt" ]
then
	echo "\tPython : Simplify parental putative_SNP list to contain only those loci with an allelic ratio on range [0.25 .. 0.75]." >> $logName;
	echo "\t\tDone." >> $logName;
	echo "\tPython : Simplify child putative_SNP list to contain only those loci with an allelic ratio on range [0.25 .. 0.75] in the parent dataset." >> $logName;
	echo "\t\tDone." >> $logName;
else
	echo "\tPython : Simplify parental putative_SNP list to contain only those loci with an allelic ratio on range [0.25 .. 0.75]." >> $logName;
	$python_exec $main_dir"scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_parent.py"            $genome $genomeUser $project $user $projectParent $projectParentUser $main_dir > $projectDirectory"trimmed_SNPs_v4.parent.txt";
	echo "\t\tDone." >> $logName;

	echo "\tPython : Simplify child putative_SNP list to contain only those loci with an allelic ratio on range [0.25 .. 0.75] in the parent dataset." >> $logName;
	$python_exec $main_dir"scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_parent_in_child.3.py" $genome $genomeUser $project $user $main_dir > $projectDirectory"trimmed_SNPs_v4.txt";
	echo "\t\tDone." >> $logName;
fi
if [ $hapmapInUse = 1 ]
then
	if [ -f $projectDirectory"trimmed_SNPs_v5.txt" ]
	then
		echo "\tPython : Simplify child putative_SNP list to contain only those loci found in the haplotype map." >> $logName;
		echo "\t\tDone." >> $logName;
	else
		echo "\tPython : Simplify child putative_SNP list to contain only those loci found in the haplotype map." >> $logName;
		$python_exec $main_dir"scripts_seqModules/putative_SNPs_from_hapmap_in_child.py" $genome $genomeUser $project $user $hapmap $hapmapUser $main_dir > $projectDirectory"trimmed_SNPs_v5.txt"
		echo "\t\tDone." >> $logName;
	fi
fi

echo "Pileup processing is complete." >> $condensedLog;
echo "\n\tPileup processing complete.\n" >> $logName;

if [ $hapmapInUse = 0 ]
then
	echo "\nPassing processing on to 'scripts_seqModules/scripts_ddRADseq/project.ddRADseq.install_4.sh' for final analysis.\n" >> $logName;
	echo   "============================================================================\n" >> $logName;
	sh $main_dir"scripts_seqModules/scripts_ddRADseq/project.ddRADseq.install_4.sh" $user $project;
else
	echo "\nPassing processing on to 'scripts_seqModules/scripts_ddRADseq/project.ddRADseq.hapmap.install_4.sh' for final analysis.\n" >> $logName;
	echo   "===================================================================================\n" >> $logName;
	sh $main_dir"scripts_seqModules/scripts_ddRADseq/project.ddRADseq.hapmap.install_4.sh" $user $project $hapmap;
fi
