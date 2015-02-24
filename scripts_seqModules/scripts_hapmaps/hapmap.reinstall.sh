#!/bin/bash -e
#
# Initialization of genome into pipeline
#   $1 : user
#   $2 : referencePloidy
#   $3 : project1 (parent or parent1)
#   $4 : project2 (child or parent2)
#   $5 : hapmap
#   $6 : main_dir

## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
referencePloidy=$2;
project1=$3;
project2=$4;
hapmap=$5;
main_dir=$(pwd)"/";

# load local installed program location variables.
. $main_dir/sh/local_installed_programs.sh;

user='darren1';
referencePloidy="2.0";
project1='SC5314_A21-s02-m09-r07';
project2='12353_A21-s02-m09-r07';
hapmap='testing';
main_dir='/heap/hapmap/bermanlab/';

echo "";
echo "Inputs to : hapmap.install_4.sh";
echo "    user            = "$user;
echo "    referencePloidy = "$referencePloidy;
echo "    project1        = "$project1;
echo "    project2        = "$project2;
echo "    hapmap          = "$hapmap;
echo "    main_dir        = "$main_dir;
echo "";
##==============================================================================
## Define locations and names to be used later.
##------------------------------------------------------------------------------
# Determine location of hapmap.
if [ -d $main_dir"users/"$user"/hapmaps/"$hapmap"/" ]
then
	hapmapDirectory=$main_dir"users/"$user"/hapmaps/"$hapmap"/";
	hapmapUser=$user;
elif [ -d $main_dir"users/default/hapmaps/"$hapmap"/" ]
then
	hapmapDirectory=$main_dir"users/default/hapmaps/"$hapmap"/";
	hapmapUser="default";
fi
# Setup process_log.txt file.
logName=$hapmapDirectory"process_log.txt";
condensedLog=$hapmapDirectory"condensed_log.txt";

echo "\thapmapDirectory = '"$hapmapDirectory"'" >> $logName;
echo "#.............................................................................." >> $logName;
echo "Running 'sh/hapmap.install_4.sh'" >> $logName;
echo "Variables passed via command-line from 'php/hapmap.install_3.php' :" >> $logName;
echo "\tuser     = '"$user"'" >> $logName;
if [ "referencePloidy" = "2" ]
then
	echo "\tparent   = '"$project1"'" >> $logName;
	echo "\tchild    = '"$project2"'" >> $logName;
else
	echo "\tparent1  = '"$project1"'" >> $logName;
	echo "\tparent2  = '"$project2"'" >> $logName;
fi
echo "\tmain_dir = '"$main_dir"'" >> $logName;
echo "" >> $logName;
echo "#=====================================#" >> $logName;
echo "# Setting up locations and variables. #" >> $logName;
echo "#=====================================#" >> $logName;
echo "Setting up for processing." >> $condensedLog;


# Determine location of project1.
if [ -d $main_dir"users/"$user"/projects/"$project1"/" ]
then
	project1Directory=$main_dir"users/"$user"/projects/"$project1"/";
	project1User=$user;
elif [ -d $main_dir"users/default/projects/"$project1"/" ]
then
	project1Directory=$main_dir"users/default/projects/"$project1"/";
	project1User="default";
fi
echo "\tproject1Directory = '"$project1Directory"'" >> $logName;

# Determine location of project2.
if [ -d $main_dir"users/"$user"/projects/"$project2"/" ]
then
	project2Directory=$main_dir"users/"$user"/projects/"$project2"/";
	project2User=$user;
elif [ -d $main_dir"users/default/projects/"$project2"/" ]
then
	projcet2Directory=$main_dir"users/default/projects/"$project2"/";
	project2User="default";
fi
echo "\tproject2Directory = '"$project2Directory"'" >> $logName;


# Get genome name from project1's "genome.txt" file.
# ...both project1 and project2 will have the same genome, as only projects matching the genome
# chosen for the hapmap are given as selection options.
genome=$(head -n 1 $project1Directory"genome.txt");
echo "\t'genome.txt' file entry." >> $logName;
echo "\t\tgenome = '"$genome"'" >> $logName;

# Determine location of project1 genome.
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

echo "\treferencePloidy = "$referencePloidy >> $logName;
if [ "$referencePloidy" = "2" ]
then
	##==============================================================================
	## For a diploid reference hapmap : Move parent SNP data to hapmap directory and
	## preprocess it for analysis.
	##------------------------------------------------------------------------------
	echo "Move SNP data files to hapmap directory." >> $logName;
	if [ ! -f $hapmapDirectory"SNPdata_parent.txt" ]
	then
		echo "\tCopy parent : 'putative_SNPs_v4.txt'" >> $logName;
		echo "\t\t to : '"$hapmapDirectory"SNPdata_parent.txt'" >> $logName;
		cp $project1Directory"putative_SNPs_v4.zip" $hapmapDirectory"SNPdata_parent.zip";
		echo "\tDecompressing parent data." >> $logName;
		cd $hapmapDirectory;
		unzip -j SNPdata_parent.zip;
		rm SNPdata_parent.zip;
		mv putative_SNPs_v4.txt SNPdata_parent.txt;
		cd $main_dir;

		# Process parent SNP file 'SNPdata_parent.txt' into condensed het SNP information.
		$python_exec $main_dir"py/hapmap.preprocess_parent.py" $genome $genomeUser $project1 $project1User $hapmap $hapmapUser $main_dir hapmap > $hapmapDirectory"SNPdata_parent.temp.txt"
		rm $hapmapDirectory"SNPdata_parent.txt"
		mv $hapmapDirectory"SNPdata_parent.temp.txt" $hapmapDirectory"SNPdata_parent.txt"
	else
		echo "\tParent data already preprocessed for use in hapmap." >> $logName;
	fi

	##==============================================================================
	## Read in 'haplotypeMap.txt' and output fragment definition files.
	##------------------------------------------------------------------------------
	haplotypeMap_file=$hapmapDirectory"haplotypeMap.txt";
	$python_exec $main_dir"py/hapmap.expand_definitions.py" $user $hapmap $main_dir

	##==============================================================================
	## Deal with installing and processing child/project2 datasets.
	##------------------------------------------------------------------------------
	# Copy child SNP dataset to hapmap directory.
	echo "\tDetermining number of child datasets in hapmap." >> $logName;
	childNum=0;
	while [ -f $hapmapDirectory"SNPdata_child."$childNum".txt" ]
	do
		childNum=`expr $childNum + 1`;
	done
	echo "\tCopy child : 'SNP_CNV_v1.zip'" >> $logName;
	echo "\t\t from : '"$project2Directory"SNP_CNV_v1.zip'" >> $logName;
	echo "\t\t to   : '"$hapmapDirectory"SNPdata_child."$childNum".zip'" >> $logName;
	cp $project2Directory"SNP_CNV_v1.zip" $hapmapDirectory"SNPdata_child."$childNum".zip";
	echo "\tDecompressing child data." >> $logName;
	cd $hapmapDirectory;
	unzip -j "SNPdata_child."$childNum".zip";
	rm "SNPdata_child."$childNum".zip";
	mv SNP_CNV_v1.txt "SNPdata_child."$childNum".txt";
	cd $main_dir

	# Process child dataset vs parental SNPs and haplotype map definitions.
	$python_exec $main_dir"py/hapmap.process_child.py" $genome $genomeUser $project2 $user $hapmap $main_dir $childNum > $hapmapDirectory"SNPdata_parent.temp.txt"

	# Delete original 'SNPdata_parent.txt' file.
	rm $hapmapDirectory"SNPdata_parent.txt";
	# Then copy 'SNPdata_parent.temp.txt' to 'SNPdata_parent.txt', as it now contains phasing information from child.
	mv $hapmapDirectory"SNPdata_parent.temp.txt" $hapmapDirectory"SNPdata_parent.txt";
	# Then delete 'SNPdata_child.0.txt' file, as no longer needed here.
	rm $hapmapDirectory"SNPdata_child."$childNum".txt";
else
	##==============================================================================
	## For a haploid reference hapmap : Move parent SNP data to hapmap directory and
	## preprocess it for analysis.
	##------------------------------------------------------------------------------
	echo "Move SNP data files to hapmap directory." >> $logName;
	if [ ! -f $hapmapDirectory"SNPdata_parent.txt" ]
	then
		echo "\tCopy parent1 : 'SNP_CNV_v1.txt'" >> $logName;
		echo "\t\t to : '"$hapmapDirectory"SNPdata_parent1.txt'" >> $logName;
		cp $project1Directory"SNP_CNV_v1.zip" $hapmapDirectory"SNPdata_parent1.zip";
		echo "\tDecompressing parent1 data." >> $logName;
		cd $hapmapDirectory;
		unzip -j SNPdata_parent1.zip;
		rm SNPdata_parent1.zip;
		mv SNP_CNV_v1.txt SNPdata_parent1.txt;
		cd $main_dir;

		echo "\tCopy parent2 : 'SNP_CNV_v1.txt'" >> $logName;
		echo "\t\t to : '"$hapmapDirectory"SNPdata_parent2.txt'" >> $logName;
		cp $project2Directory"SNP_CNV_v1.zip" $hapmapDirectory"SNPdata_parent2.zip";
		echo "\tDecompressing parent2 data." >> $logName;
		cd $hapmapDirectory;
		unzip -j SNPdata_parent2.zip;
		rm SNPdata_parent2.zip;
		mv SNP_CNV_v1.txt SNPdata_parent2.txt;
		cd $main_dir;

		# Process parent SNP files 'SNPdata_parent1.txt' & 'SNPdata_parent2.txt' into condensed het SNP information.
		$python_exec $main_dir"py/hapmap.preprocess_haploid_parents.py" $genome $genomeUser $project1 $project1User $project2 $project2User $hapmap $hapmapUser $main_dir hapmap > $hapmapDirectory"SNPdata_parent.temp.txt"
		mv $hapmapDirectory"SNPdata_parent.temp.txt" $hapmapDirectory"SNPdata_parent.txt"

		# Delete original 'SNPdata_parent1.txt' and 'SNPdata_parent2.txt' file, no longer needed here.
		rm $hapmapDirectory"SNPdata_parent1.txt";
		rm $hapmapDirectory"SNPdata_parent2.txt";
	else
		echo "\tParent1 & parent2 data already preprocessed for use as hapmap." >> $logName;
	fi

	## Generate "complete.txt" to indicate processing has completed normally.
	completeFile=$hapmapDirectory"complete.txt";
	timestamp=$(date +%T);
	echo $timestamp > $completeFile;
	echo "\tGenerated 'complete.txt' file." >> $logName;
	chmod 0755 $completeFile;
fi

echo "Concluding analysis." >> $condensedLog;

## Delete 'working.txt' file to let pipeline know that processing has completed, but hapmap is available for additional entries.
rm $main_dir"users/"$user"/hapmaps/"$hapmap"/working.txt";


##==============================================================================
## Cleanup intermediate files.
##------------------------------------------------------------------------------
rm $hapmapDirectory"condensed_log.txt";
rm $hapmapDirectory"process_log.txt";
