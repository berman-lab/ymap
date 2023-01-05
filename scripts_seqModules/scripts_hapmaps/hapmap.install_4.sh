#!/bin/bash -e

## All created files will have permission 760
umask 007;

### define script file locations.
user=$1;
referencePloidy=$2;
project1=$3;
project2=$4;
hapmap=$5;

main_dir=$(pwd)"/../../";

# load local installed program location variables.
. $main_dir/local_installed_programs.sh;


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
logName=$hapmapDirectory"process_log.txt";
condensedLog=$hapmapDirectory"condensed_log.txt";
echo "" >> $logName;
echo "Running 'scripts_seqModules/scripts_hapmaps/hapmap.install_4.sh'" >> $logName;
echo "Variables passed via command-line from 'scripts_seqModules/scripts_hapmaps/hapmap.install_3.php' :" >> $logName;
echo "    user                        = "$user >> $logName;
echo "    referencePloidy             = "$referencePloidy >> $logName;
echo "    project1 (parent)           = "$project1 >> $logName;
echo "    project2 (child)            = "$project2 >> $logName;
echo "    hapmap                      = "$hapmap >> $logName;
echo "    main_dir                    = "$main_dir >> $logName;
echo "#.............................................................................." >> $logName;
echo "" >> $logName;
echo "#=====================================#" >> $logName;
echo "# Setting up locations and variables. #" >> $logName;
echo "#=====================================#" >> $logName;
echo "Setting up for processing." >> $condensedLog;
echo "Important variables :" >> $logName;
echo "    hapmap user                 = '"$hapmapUser"'" >> $logName;
echo "    hapmap directory            = '"$hapmapDirectory"'" >> $logName;

# Determine location of project1 (parent).  Is it in user or default account?
if [ -d $main_dir"users/"$user"/projects/"$project1"/" ]
then
	project1Directory=$main_dir"users/"$user"/projects/"$project1"/";
	project1User=$user;
elif [ -d $main_dir"users/default/projects/"$project1"/" ]
then
	project1Directory=$main_dir"users/default/projects/"$project1"/";
	project1User="default";
fi
echo "    project1 (parent) directory = '"$project1Directory"'" >> $logName;

# Determine location of project2 (child).  Is it in user or default account?
if [ -d $main_dir"users/"$user"/projects/"$project2"/" ]
then
	project2Directory=$main_dir"users/"$user"/projects/"$project2"/";
	project2User=$user;
elif [ -d $main_dir"users/default/projects/"$project2"/" ]
then
	projcet2Directory=$main_dir"users/default/projects/"$project2"/";
	project2User="default";
fi
echo "    project2 (child) directory  = '"$project2Directory"'" >> $logName;

# Get genome name from project1's "genome.txt" file.
# ...both project1 and project2 will have the same genome, as only projects matching the genome
# chosen for the hapmap are given as selection options.
genome=$(head -n 1 $project1Directory"genome.txt");
echo "    genome                      = '"$genome"'" >> $logName;

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
echo "    genome directory            = '"$genomeDirectory"'" >> $logName;

# Get reference FASTA file name from "reference.txt";
genomeFASTA=$(head -n 1 $genomeDirectory"reference.txt");
echo "    genome FASTA file           = '"$genomeFASTA"'" >> $logName;
echo "    reference ploidy            = '"$referencePloidy"'" >> $logName;

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
		$python_exec $main_dir"scripts_seqModules/scripts_hapmaps/hapmap.preprocess_parent.py" $genome $genomeUser $project1 $project1User $hapmap $hapmapUser $main_dir hapmap > $hapmapDirectory"SNPdata_parent.temp.txt"
		rm $hapmapDirectory"SNPdata_parent.txt"
		mv $hapmapDirectory"SNPdata_parent.temp.txt" $hapmapDirectory"SNPdata_parent.txt"
	else
		echo "\tParent data already preprocessed for use in hapmap." >> $logName;
	fi

	##==============================================================================
	## Read in 'haplotypeMap.txt' and output fragment definition files.
	##------------------------------------------------------------------------------
	$python_exec $main_dir"scripts_seqModules/scripts_hapmaps/hapmap.expand_definitions.py" $user $hapmap $main_dir

	##==============================================================================
	## Deal with installing and processing child/project2 datasets.
	##------------------------------------------------------------------------------
	# Copy child SNP dataset to hapmap directory.
	echo "\tDetermining number of child datasets in hapmap." >> $logName;
	childNum=0;
	while [ -f $hapmapDirectory"haplotypeFragments."$childNum".txt" ]
	do
		childNum=`expr $childNum + 1`;
	done
	echo "\t\tThere are "$childNum" child datasets used in hapmap." >> $logName;
	childNum=`expr $childNum - 1`;   # at least one will always be found, but counting of map entries is zero based.
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
	$python_exec $main_dir"scripts_seqModules/scripts_hapmaps/hapmap.process_child.py" $genome $genomeUser $project2 $user $hapmap $main_dir $childNum > $hapmapDirectory"SNPdata_parent.temp.txt"

	# Delete original 'SNPdata_parent.txt' file.
	rm $hapmapDirectory"SNPdata_parent.txt";
	# Then copy 'SNPdata_parent.temp.txt' to 'SNPdata_parent.txt', as it now contains phasing information from child.
	mv $hapmapDirectory"SNPdata_parent.temp.txt" $hapmapDirectory"SNPdata_parent.txt";
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
		$python_exec $main_dir"scripts_seqModules/scripts_hapmaps/hapmap.preprocess_haploid_parents.py" $genome $genomeUser $project1 $project1User $project2 $project2User $hapmap $hapmapUser $main_dir hapmap > $hapmapDirectory"SNPdata_parent.temp.txt"
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
