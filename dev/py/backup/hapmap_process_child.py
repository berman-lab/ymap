# Input Arguments
#   1) user     : 'darren'
#   2) project  : 'test_Ca'
#   3) genome   : 'Candida_albicans_SC5314_etc'
#   4) main_dir : '/home/bermanj/shared/links/'
#	5) logName  :
#
# Process input files:
#	1) Raw SNP data                 : $workingDir"users/"$user"/projects/"$project"/putative_SNPs_v4.txt".
#	2) FASTA file name              : $workingDir"users/default/genomes/default/reference.txt",
#	                               or $workingDir"users/"$user"/genomes/default/reference.txt" as $FastaName.
#   3) hapmapParent_preprocessed_SNPs.txt
#           [P_chr_name][P_position][Allele1][Alleles2]

# Generate output file:
#	1) a simplified pileup file containing only heterozygous loci in parental dataset.
#		0) chr_name      : Chromosome name string.
#		1) bp_coordinate : bp coordinate along chromosome.
#		2) baseCall_1    : alphabetical first call of heterozgous locus.
#		3) baseCall_2    : alphabetical second call of heterozygous locus.
#		Comment lines in output begin with '#'.
#

import string, sys, re, time, os
genome            = sys.argv[1]
genomeUser        = sys.argv[2]
projectParent     = sys.argv[3]
projectParentUser = sys.argv[4]
main_dir          = sys.argv[7]
logName           = sys.argv[8]
inputFile1        = main_dir+"users/"+projectParentUser+"/projects/"+projectParent+"/putative_SNPs_v4.txt"

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\n*--------------------------------------------------*")
	myfile.write("\n| Log of 'py/hapmap_process_parent.py'             |")
	myfile.write("\n*--------------------------------------------------*\n")

#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+genomeUser+"/genomes/"+genome+"/"

with open(logName, "a") as myfile:
	myfile.write("\n\tProcessing standard bin fragmented genome file.")

#============================================================================================================
# Load FastaName from 'reference.txt' for genome in use.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\tIdentifying name of reference FASTA file.\n")
reference_file = genomeDirectory + '/reference.txt'
refFile        = open(reference_file,'r')
FastaName      = refFile.read().strip()
refFile.close()
FastaName      = FastaName.replace(".fasta", "")

#============================================================================================================
# Process 'hapmapParent_preprocessed_SNPs.txt' file for projectParent to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\tProcessing datasetParent 'preprocessed_SNPs.txt' file -> het loci.\n")

# Look up chromosome name strings for genome in use.
#     Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/figure_definitions.txt"
figureDefinition_file  = genomeDirectory + 'figure_definitions.txt'
figureDefinitionFile   = open(figureDefinition_file,'r')
figureDefinitionData   = figureDefinitionFile.readlines()

# Example lines in figureDefinition_file:
#     Chr  Use   Label   Name                         posX   posY   width   height
#     1    1     Chr1    Ca21chr1_C_albicans_SC5314   0.15   0.8    0.8     0.0625
#     2    1     Chr2    Ca21chr2_C_albicans_SC5314   0.15   0.7    *       0.0625
#     0    0     Mito    Ca19-mtDNA                   0.0    0.0    0.0     0.0
with open(logName, "a") as myfile:
	myfile.write("\tDetermining number of chromosomes of interest in genome.\n")

# Determine the number of chromosomes of interest in genome.
chrName_maxcount = 0
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line))
	chr_num = line_parts[0]
	if chr_num.isdigit():
		chr_num    = int(float(line_parts[0]))
		chr_use    = int(float(line_parts[1]))
		chr_label  = line_parts[2]
		chr_name   = line_parts[3]
		if chr_num > 0:
			if chr_num > chrName_maxcount:
				chrName_maxcount = chr_num
figureDefinitionFile.close()

# Pre-allocate chrName_array
chrName = []
for x in range(0, chrName_maxcount):
	chrName.append([])

with open(logName, "a") as myfile:
	myfile.write("\n\tGathering name strings for chromosomes.")

# Gather name strings for chromosomes, in order.
figureDefinitionFile  = open(figureDefinition_file,'r')
chrCounter = 0;
chrNums    = [];
chrNames   = [];
chrLabels  = [];
chrShorts  = [];
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line))
	chr_num = line_parts[0]
	if chr_num.isdigit():
		chr_num                        = int(float(line_parts[0]))
		chrNums.append(chr_num);
		chrCounter += chrCounter;
		chr_use                        = int(float(line_parts[1]))
		chr_label                      = line_parts[2]
		chrLabels.append(chr_label);
		chr_name                       = line_parts[3]
		chrNames.append(chr_name);
		chr_nameShort                  = chr_label
		chrShorts.append(chr_nameShort);
		if chr_num != 0:
			chrName[int(float(chr_num))-1] = chr_name
			with open(logName, "a") as myfile:
				myfile.write("\n\t\t" + str(chr_num) + " : " + chr_name + " = " + chr_nameShort)
figureDefinitionFile.close()

# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount
with open(logName, "a") as myfile:
	myfile.write("\n\t\tMax chr string : "+str(chrCount)+"\n")
#............................................................................................................

with open(logName, "a") as myfile:
	myfile.write("\n\tOpen datasetParent 'putative_SNPs_v4.txt' file.\n")

# Open dataset 'putative_CNVs_v1.txt' file.
print '### InputFile1 = ' + inputFile1
data          = open(inputFile1,'r')
#............................................................................................................

count            = 0
old_chr          = 0
fragment_found   = 0
last_fragment    = 0
current_fragment = 0
log_offset       = 0

print '### Number of Chromosomes = ' + str(chrCount)
for x in range(0,chrCount):
	if (chrNums[x] != 0):
		print '### \t' + str(x+1) + ' : ' + str(chrName[x])

with open(logName, "a") as myfile:
	myfile.write("\tGathering read coverage data for each fragment.")

# Process parent 'putative_SNPs_v4.txt' file, line by line.
for line in data:
	# example lines from SNP pileup file:
	#		chromosome                   coord   ref   A    T    G    C
	#		Ca21chr1_C_albicans_SC5314   2461    T     0    108  1    0
	#		Ca21chr1_C_albicans_SC5314   2463    A     110  1    0    0
	#		Ca21chr1_C_albicans_SC5314   2464    C     0    0    1    108
	count += 1
	parentLine = string.strip(line)
	parentLine = parentLine.split('\t')
	P_chr_name = parentLine[0]   # chr name of bp.			: Ca21chrR_C_albicans_SC5314
	P_position = parentLine[1]   # chr position of bp.		: 2286371
	P_refBase  = parentLine[2]   # reference base at bp.	: T
	P_countA   = parentLine[3]   # count of A.				: 100
	P_countT   = parentLine[4]   # count of T.				: 0
	P_countG   = parentLine[5]   # count of G.				: 0
	P_countC   = parentLine[6]   # count of C.				: 1

	# Determine if parent data point is heterozygous.
	P_counts   = [int(float(P_countA)), int(float(P_countT)), int(float(P_countG)), int(float(P_countC))]
	maxAllele1 = max(P_counts)
	P_het      = False
	P_ratio    = maxAllele1/float(sum(P_counts))
	if (P_ratio > 0.25) and (P_ratio < 0.75):
		P_het = True
	else:
		P_het = False

	if P_het:	# ...then output a pileup line.
		# Identify which chromosome this data point corresponds to.
		chr = 0
		for x in range(0,chrCount):
			if (chrNums[x] != 0):
				if chrName[x] == P_chr_name:
					chr = x+1

		# Identify base calls which contribute to heterozygous call.
		alleleData    = [(int(float(P_countA), 'A'), (int(float(P_countT), 'T'), (int(float(P_countG), 'G'), (int(float(P_countC), 'C')]
		sortedAlleles = sorted(alleleData, key=lambda alleleDatum: alleleDatum[0])
		Alleles       = [sortedAlleles[3], sortedAlleles[2]];
		Alleles.sort()

		# Display status updates to log file.
		if old_chr != chr:
			print '### chr change : ' + str(old_chr) + ' -> ' + str(chr)
			with open(logName, "a") as myfile:
				myfile.write("\n\t    " + str(old_chr) + " -> " + str(chr) + " = " + P_chr_name + "\n")
				myfile.write("\t\t1........01........01........01........01........01........01........01........01........01........0")
		if chr!=0:
			# Reset for each new chromosome examined.
			if old_chr != chr:
				if log_offset != 0:
					log_offset_string = " "*((log_offset)%100)
					with open(logName, "a") as myfile:
						myfile.write("\n\t\t" + log_offset_string)
				count = 1
			else:
				count += 1

			# display status updates to log file.
			if (count%50000 == 0):
				log_offset += 1
				if (count%50000*50) == 0:
					if (current_fragment-1) == 0:
						with open(logName, "a") as myfile:
							myfile.write("\n\t\t")
					else:
						with open(logName, "a") as myfile:
							myfile.write(" " + str(current_fragment-1) + "\n\t\t")
				with open(logName, "a") as myfile:
					myfile.write(".")

		# Output new pileup line
		#	0) chr_name      : Chromosome name string.
		#	1) bp_coordinate : bp coordinate along chromosome.
		#	2) baseCall_1    : alphabetical first call of heterozgous locus.
		#	3) baseCall_2    : alphabetical second call of heterozygous locus.
		print P_chr_name + '\t' + str(P_position) + '\t' + Alleles[0] + '\t' + Alleles[1]


		# Reset old_chr/last_fragment to current coordinate chromosome before moving to next line in pileup. 
		old_chr       = chr

print "### ", time.clock() - t1, "seconds to parse project SNP data."
t2 = time.clock()
print '### Data from each fragment: [chrNum, bpStart, bpEnd, Max, Ave, Length]'

#------------------------------------------------------------------------------------------------------------
# End of code section to output information about parental heterozygous loci. 
#============================================================================================================

print "### ", time.clock() - t1, "seconds to output basic stats of each restriction fragment."
print "### ", time.clock() - t0, "seconds to complete processing of fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\n\tTime to process = " + str(time.clock()-t0) )
	myfile.write("\n\t* 'py/hapmap_process_parent.py' completed. *\n")
