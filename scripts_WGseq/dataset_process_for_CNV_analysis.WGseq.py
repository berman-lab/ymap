# Input Arguments
#   1) user       : 'darren'
#   2) project    : 'test_Ca'
#   3) genome     : 'Candida_albicans_SC5314_etc'
#   4) genomeUser : ''
#   5) main_dir   : '/home/bermanj/shared/links/'
#	6) logName    :
#
# Process input files:
#	1) Raw CNV data                 : $workingDir"users/"$user"/projects/"$project"/SNP_CNV_v1.txt".
#	2) FASTA file name              : $workingDir"users/default/genomes/default/reference.txt",
#	                               or $workingDir"users/"$user"/genomes/default/reference.txt" as $FastaName.
#	3) Coordinates of standard bins : $workingDir"users/default/genomes/"$genome"/"$FastaName".standard_bins.fasta",
#	                               or $workingDir"users/"$user"/genomes/"$genome"/"$FastaName".standard_bins.fasta".

# Generate output file:
#	1) a simplified pileup file containing average read counts per standard bin.   [chr_num,bp_start,bp_end, data_ave]
#		0) chr_num   : Numerical chromosome identifier, defined for each genome in "figure_details.txt".
#		1) bp_start  : Start bp coordinate along chromosome.
#		2) bp_end    : End bp coordinate along chromosome.
#		3) reads_ave : Average reads across fragment.
#		Comment lines in output begin with '#'.
#

import string, sys, re, time
user       = sys.argv[1]
project    = sys.argv[2]
genome     = sys.argv[3]
genomeUser = sys.argv[4]
main_dir   = sys.argv[5]
logName    = sys.argv[6]
inputFile  = main_dir+"users/"+user+"/projects/"+project+"/SNP_CNV_v1.txt"

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\t\t\t*======================================================================*\n")
	myfile.write("\t\t\t| Log of 'dataset_process_for_CNV_analysis.WGseq.py'                   |\n")
	myfile.write("\t\t\t*----------------------------------------------------------------------*\n")

#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+genomeUser+"/genomes/"+genome+"/"

with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\n")
	myfile.write("\t\t\t|\tProcessing standard bin fragmented genome file.\n")

#============================================================================================================
# Load FastaName from 'reference.txt' for genome in use.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tIdentifying name of reference FASTA file.\n")
reference_file = genomeDirectory + '/reference.txt'
refFile        = open(reference_file,'r')
FastaName      = refFile.read().strip()
refFile.close()
FastaName      = FastaName.replace(".fasta", "")

#============================================================================================================
# Process standard-bin fragmented reference genome file.
#------------------------------------------------------------------------------------------------------------
# Example fragmented FASTA header line.
#     >Ca_a.chr1 (9638..10115) (478bp) [*]
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tLoading standard bin fragmented genome FASTA file.\n")

# Open restriction-digested genome FASTA file.
standardBins_FASTA_file = genomeDirectory + FastaName + ".standard_bins.fasta"
standardBins_FASTA_data = open(standardBins_FASTA_file,'r')

# Setup array and counter for tracking fragment definition data.
fragments        = []
fragment_counter = 0

## Process digested FASTA genome file, line by line.
while True:
	# Line pairs have the following structure.
	#    >Ca_a.chr1 (9638..10115) (478bp) [*]
	#    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
	line1 = standardBins_FASTA_data.readline()
	line2 = standardBins_FASTA_data.readline()
	if not line2:
		break  # EOF
	first_char = line1[:1]
	if first_char == ">":
		# Line is header to FASTA entry.
		line_parts             = string.split(string.strip(line1))
		chrGenomeAndNum_string = line_parts[0]
		bp_coordinate_string   = line_parts[1]
		fragment_size_string   = line_parts[2]
		if len(line_parts) > 3:
			fragment_usable_string = line_parts[3]
			if fragment_usable_string[1] == "*":
				# Fragment is usable, so the details should be placed into fragments structure.
				# split the chr string by '.' character, then trim off the first three characters ('chr') from the second substring.
				#   string has format of : ">Ca_a.chr1"
				genomeName_string,chrNum_string = chrGenomeAndNum_string.split(".")
				chr_num                         = int(float(chrNum_string.replace("chr","")))

				#   string has format of : "(9638..10115)"
				coordinates = bp_coordinate_string.replace('(','').replace(')','').replace('..',' ').split()
				bp_start    = int(float(coordinates[0]))
				bp_end      = int(float(coordinates[1]))
				reads_count = 0   # placeholder value.
				reads_max   = 0   # placeholder value.
				reads_ave   = 0   # placeholder value.

				fragments.append([chr_num,bp_start,bp_end,reads_count,reads_max,reads_ave])
				fragment_counter += 1
standardBins_FASTA_data.close()

with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tThe standard bin fragmented genome FASTA file has been loaded.\n\t\t\t|")

# Put fragment counter into a general use variable.
numFragments = fragment_counter
#------------------------------------------------------------------------------------------------------------
# End of code section to parse restriction fragments from genome.
#============================================================================================================

print "### ", time.clock() - t0, "seconds to parse restriction fragments from digested genome."
t1 = time.clock()
print "### Starting read count data processing."

#============================================================================================================
# Process 'SNP_CNV_v1.txt' file to determine read count max and average.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t|\tProcessing dataset 'SNP_CNV_v1.txt' file -> max and average read counts per fragment.\n")

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
	myfile.write("\t\t\t|\tDetermining number of chromosomes of interest in genome.\n")

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
	myfile.write("\t\t\t|\tGathering name strings for chromosomes.\n")

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
				myfile.write("\t\t\t|\t" + str(chr_num) + " : " + chr_name + " = " + chr_nameShort + "\n")
figureDefinitionFile.close()

# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tMax chr string : "+str(chrCount)+"\n")
#............................................................................................................

with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tOpen dataset 'SNP_CNV_v1.txt' file.\n")

# Open dataset 'SNP_CNV_v1.txt' file.
print '### InputFile = ' + inputFile
datafile      = inputFile;
data          = open(datafile,'r')
#............................................................................................................

count            = 0
old_chr          = 0
fragment_found   = 0
last_fragment    = 0
current_fragment = 0
log_count        = 0
log_offset       = 0

print '### Number of Chromosomes = ' + str(chrCount)
for x in range(0,chrCount):
	if (chrNums[x] != 0):
		print '### \t' + str(x+1) + ' : ' + str(chrName[x])

print "###" + str(numFragments)

with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tGathering read coverage data for each fragment.\n\t\t\t|")

# Process 'SNP_CNV_v1.txt' file, line by line.
for line in data:
	# example lines from CNV pileup file:
	#     chromosomeNam                 pos             totalReads        0
	#     Ca21chr1_C_albicans_SC5314    2388924         123               0
	#     Ca21chr1_C_albicans_SC5314    2388925         135               0
	count += 1
	line_parts = string.split(string.strip(line))
	chr_name   = line_parts[0]   # chr name of bp.		: Ca21chrR_C_albicans_SC5314
	position   = line_parts[1]   # chr position of bp.	: 2286371
	readCount  = line_parts[2]   # read count at bp.	: 12

	# Attempt to match up current data line with pre-determined restriction fragments.
	found = 0

	# Identify which chromosome this data point corresponds to.
	chr = 0
	for x in range(0,chrCount):
		#if (chrNums[x] != 0):
			#	print str(chrName[x])
		if (chrNums[x] != 0):
			if chrName[x] == chr_name:
				chr = x+1

	# Convert to integers.
	pos_point  = int(position)
	data_point = float(readCount)

	if old_chr != chr:
		print '### chr change : ' + str(old_chr) + ' -> ' + str(chr)
		with open(logName, "a") as myfile:
			myfile.write("\n\t\t\t|\n\t\t\t|\t" + str(old_chr) + " -> " + str(chr) + " = " + chr_name + "\n")
			myfile.write("\t\t\t|\t1........01........01........01........01........01........01........01........01........01........0\n")

	if chr!=0:
		# Reset for each new chromosome examined.
		if old_chr != chr:
			if log_offset != 0:
				log_offset_string = " "*((log_offset)%100)
				with open(logName, "a") as myfile:
					myfile.write("\t\t\t|\t" + log_offset_string)
			count            = 1
			fragment_found   = 0
			current_fragment = 0

		# If (fragment_found == 0), look to see if current coordinate matches any defined fragments.
		if fragment_found == 0:
			for frag in range(current_fragment,numFragments):
				# Check if current coordinate is consistent with this fragment : fragments[frag-1] = [chr_num,bp_start,bp_end,readSum,readMax,readAve]
				# print str(chr)+": "+str(pos_point)+"; "+str(fragments[frag-1][0])+":"+str(fragments[frag-1][1])+"..."+str(fragments[frag-1][2])
				if chr == fragments[frag-1][0] and pos_point >= fragments[frag-1][1] and pos_point <= fragments[frag-1][2]:
					fragment_found   = 1
					current_fragment = frag
					break

		# If (fragment_found == 1), add current bp coordinate data to fragment data.
		if fragment_found == 1:
			if last_fragment != current_fragment:
				log_count  += 1
				log_offset += 1
				if ((log_count-1)%100) == 0:
					if (log_count-1) == 0:
						with open(logName, "a") as myfile:
							myfile.write("\t\t\t|\t")
					else:
						with open(logName, "a") as myfile:
							myfile.write(" " + str(log_count-1) + "\n\t\t\t|\t")
				with open(logName, "a") as myfile:
					myfile.write(".")

			# Adds current coordinate repetitiveness score to fragment total count : fragments[frag-1] = [chr_num,bp_start,bp_end,readSum,readMax,readAve]
			fragments[current_fragment-1][3] += data_point

			# If current coordinate read count is highest so far for fragment, update max : fragments[frag-1] = [chr_num,bp_start,bp_end,readSum,readMax,readAve]
			if data_point > fragments[current_fragment-1][4]:
				fragments[current_fragment-1][4] = data_point

			# If current coordinate is at (or after) the end of a fragment, update fragment_found to 0 : fragments[frag-1] = [chr_num,bp_start,bp_end,readSum,readMax,readAve]
			if pos_point >= fragments[current_fragment-1][2]:
				fragment_found = 0

	# Reset old_chr to current coordinate chromosome before moving to next line in pileup. 
	old_chr       = chr
	last_fragment = current_fragment

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t|\tCalculating average read coverage per fragment.\n")

# Calculate average read coverage per each fragment.
for fragment in range(1,numFragments):
	# Calculate average read coverage for this fragment.
	#       fragments[fragment-1] = [chr_num,bp_start,bp_end, readSum,readMax,readAve]
	fragments[fragment-1][5] = fragments[fragment-1][3]/float(fragments[fragment-1][2]-fragments[fragment-1][1]+1)
#------------------------------------------------------------------------------------------------------------
# End of code section to parse read count max/average data.
#============================================================================================================


print "### ", time.clock() - t1, "seconds to process the pileup file."
t2 = time.clock()
print '### Number of fragments = ' + str(numFragments)
print '### Data from each fragment: [chrNum, bpStart, bpEnd, aveDepth]'

#============================================================================================================
# Code section to output information about read average per standard bin genome fragment.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tOutput condensed CNV per fragment information.\n")
for fragment in range(1,numFragments):
	# Output a line for each fragment.
	#     fragments[fragment-1] = [chr_num,bp_start,bp_end, aveDepth]
	#     0) chr_num
	#     1) bp_start
	#     2) bp_end
	#     3) average reads
	chr_num         = fragments[fragment-1][0]
	bp_start        = fragments[fragment-1][1]
	bp_end          = fragments[fragment-1][2]
	read_average    = fragments[fragment-1][5]
	print str(chr_num) + '\t' + str(bp_start) + '\t' + str(bp_end) + '\t' + str(read_average)
#------------------------------------------------------------------------------------------------------------
# End of code section to output information about fragments. 
#============================================================================================================


print "### ", time.clock() - t1, "seconds to output basic stats of each restriction fragment."
print "### ", time.clock() - t0, "seconds to complete processing of fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tTime to process = " + str(time.clock()-t0) +"\n")
	myfile.write("\t\t\t|\t'scripts_WGseq/dataset_process_for_CNV_analysis.WGseq.py' completed.\n")
	myfile.write("\t\t\t*----------------------------------------------------------------------*\n")
	myfile.write("\t\t\t| End of Log from 'dataset_process_for_CNV_analysis.WGseq,py'          |\n")
	myfile.write("\t\t\t*======================================================================*\n")
