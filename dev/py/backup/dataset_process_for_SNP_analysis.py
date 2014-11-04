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
#	3) Coordinates of standard bins : $workingDir"users/default/genomes/"$genome"/"$FastaName".standard_bins.fasta",
#	                               or $workingDir"users/"$user"/genomes/"$genome"/"$FastaName".standard_bins.fasta".

# Generate output file:
#	1) a simplified pileup file containing number of parental het loci that are [HOM, HET, oddHET] in child per standard bin.
#		0) chr_num  : Numerical chromosome identifier, defined for each genome in "figure_details.txt".
#		1) bp_start : Start bp coordinate along chromosome.
#		2) bp_end   : End bp coordinate along chromosome.
#		3) het from parent -> het in child.
#		4) het from parent -> oddhet in child.
#		5) het from parent -> hom in child.
#		Comment lines in output begin with '#'.
#

import string, sys, re, time
genome            = sys.argv[1]
genomeUser        = sys.argv[2]
projectParent     = sys.argv[3]
projectParentUser = sys.argv[4]
projectChild      = sys.argv[5]
projectChildUser  = sys.argv[6]
main_dir          = sys.argv[7]
logName           = sys.argv[8]
inputFile1        = main_dir+"users/"+projectParentUser+"/projects/"+projectParent+"/putative_SNPs_v4.txt"
inputFile2        = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/putative_SNPs_v4.txt"

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\n*--------------------------------------------------*")
	myfile.write("\n| Log of 'dataset_process_for_SNP_analysis.py'     |")
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
# Process standard-bin fragmented reference genome file.
#------------------------------------------------------------------------------------------------------------
# Example fragmented FASTA header line.
#     >Ca_a.chr1 (9638..10115) (478bp) [*]
with open(logName, "a") as myfile:
	myfile.write("\tLoading standard bin fragmented genome FASTA file.\n")

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
	first_char = line1[:1];
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
				coordinates  = bp_coordinate_string.replace('(','').replace(')','').replace('..',' ').split()
				bp_start     = int(float(coordinates[0]))
				bp_end       = int(float(coordinates[1]))
				HOM_count    = 0   # placeholder value.
				HET_count    = 0   # placeholder value.
				oddHET_count = 0   # placeholder value.

				fragments.append([chr_num,bp_start,bp_end,HOM_count,HET_count,oddHET_count])
				fragment_counter += 1
standardBins_FASTA_data.close()

with open(logName, "a") as myfile:
	myfile.write("\tThe standard bin fragmented genome FASTA file has been loaded.\n")

# Put fragment counter into a general use variable.
numFragments = fragment_counter
#------------------------------------------------------------------------------------------------------------
# End of code section to parse fragments from genome.
#============================================================================================================

print "### ", time.clock() - t0, "seconds to parse restriction fragments from digested genome."
t1 = time.clock()
print "### Starting read count data processing."

#============================================================================================================
# Process 'putative_SNPs_v4.txt' file for projectParent to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\tProcessing datasetParent 'putative_SNPs_v4.txt' file -> het loci.\n")

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

print "###" + str(numFragments)

with open(logName, "a") as myfile:
	myfile.write("\tGathering read coverage data for each fragment.")

# Process parent 'putative_SNPs_v4.txt' file, line by line.
for line in data:
	# example lines from CNV pileup file:
	#		chromosome                   coord   ref   A    T    G    C
	#		Ca21chr1_C_albicans_SC5314   2461    T     0    108  1    0
	#		Ca21chr1_C_albicans_SC5314   2463    A     110  1    0    0
	#		Ca21chr1_C_albicans_SC5314   2464    C     0    0    1    108
	count += 1
	parentLine   = string.strip(line)
	parentLine   = parentLine.split('\t')
	P_chr_name   = parentLine[0]   # chr name of bp.			: Ca21chrR_C_albicans_SC5314
	P_position   = parentLine[1]   # chr position of bp.		: 2286371
	P_refBase    = parentLine[2]   # reference base at bp.		: T
	P_countA     = parentLine[3]   # count of A.				: 100
	P_countT     = parentLine[4]   # count of T.				: 0
	P_countG     = parentLine[5]   # count of G.				: 0
	P_countC     = parentLine[6]   # count of C.				: 1

	# Determine if parent data point is heterozygous.
	P_counts       = [int(float(P_countA)), int(float(P_countT)), int(float(P_countG)), int(float(P_countC))]
	maxAlleleCount = max(P_counts)
	P_het          = False
	P_ratio        = maxAlleleCount/float(sum(P_counts))
	if (P_ratio > 0.25) and (P_ratio < 0.75):
		P_het = True
	else:
		P_het = False

	if P_het:
		if projectChild == projectParent:
			C_het    = True
			C_oddHet = False
		else:
			# Attempt to find a matching data line in projectChild 'putative_SNPs_v4.txt' file, using grep.
			# grep -P "Ca21chr1_C_albicans_SC5314\t2461\t" putative_SNPs_v4.txt 
			searchString = P_chr_name+"\t"+P_position+"\t"
			result = ""

			# Functional, but slow.
			for lineChild in open(inputFile2,'r'):
				if searchString in lineChild:
					result = lineChild
					break
			#	try:
			#		childLine = subprocess.check_output('grep %s '+inputFile2 % (searchString), shell = True).strip.split('\n')
			#	except subprocess.CalledProcessError as e:
			#		if e.returncode > 1:
			#			raise
			#		childLine = []
			if result=="":
				# locus is not found in child dataset => homozygous.
				C_het     = False
				#C_ratio   = "inf"
				#outString = "HOM"
			else:
				# locus is found in child dataset => heterozygous or odd het.
				C_het        = True
				childLine    = string.strip(result)
				childLine    = childLine.split('\t')
				C_chr_name   = childLine[0]
				C_position   = childLine[1]
				C_refBase    = childLine[2]
				C_countA     = childLine[3]
				C_countT     = childLine[4]
				C_countG     = childLine[5]
				C_countC     = childLine[6]

				# Determine if child data point is heterozygous.
				C_counts       = [int(float(C_countA)), int(float(C_countT)), int(float(C_countG)), int(float(C_countC))]
				maxAlleleCount = max(C_counts)
				C_ratio        = maxAlleleCount/float(sum(C_counts))
				if (C_ratio > 0.25) and (C_ratio < 0.75):
					C_oddHet  = False
					#outString = "HET"
				else:
					C_oddHet  = True
					#outString = "oddHET"

			#print "\tP|" + string.strip(line)   + "\t|" + str(P_ratio)
			#print "\tC|" + string.strip(result) + "\t|" + str(C_ratio) + "\t" + outString + "\n"

		#=====================================
		# Find fragment bin which holds data.
		#-------------------------------------
		found = 0

		# Identify which chromosome this data point corresponds to.
		chr = 0
		for x in range(0,chrCount):
			#if (chrNums[x] != 0):
			#	print str(chrName[x])
			if (chrNums[x] != 0):
				if chrName[x] == P_chr_name:
					chr = x+1

		# Convert to integers.
		pos_point  = int(float(P_position))

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
				count            = 1
				fragment_found   = 0
				current_fragment = 0

			# If (fragment_found == 0), look to see if current coordinate matches any defined fragments.
			if fragment_found == 0:
				for frag in range(current_fragment,numFragments):
					# Check if current coordinate is consistent with this fragment : fragments[frag-1] = [chr_num,bp_start,bp_end, data_count,data_max,ave_read_count, repetDataCount,repetMax,repetAve]
					if chr == fragments[frag-1][0] and pos_point >= fragments[frag-1][1] and pos_point <= fragments[frag-1][2]:
						fragment_found   = 1
						current_fragment = frag
						break

			# If (fragment_found == 1), add current bp coordinate data to fragment data.
			if fragment_found == 1:
				# display status updates to log file.
				if last_fragment != current_fragment:
					log_offset += 1
					if ((current_fragment-1)%100) == 0:
						if (current_fragment-1) == 0:
							with open(logName, "a") as myfile:
								myfile.write("\n\t\t")
						else:
							with open(logName, "a") as myfile:
								myfile.write(" " + str(current_fragment-1) + "\n\t\t")
					with open(logName, "a") as myfile:
						myfile.write(".")

				# If current coordinate is after end of a fragment, update fragment_found to 0.
				if pos_point > fragments[current_fragment-1][2]:
					fragment_found = 0
				else:
					#==========================================================
					# Add HET/oddHET/HOM status of locus to fragment data list.
					#----------------------------------------------------------
					#     C_het = False                   => homozygous.
					#     C_het = True, C_oddHet = False  => heterozygous.
					#     C_het = True, C_oddHet = True   => odd heterozygous.
					# fragment = [0=chr_num,1=bp_start,2=bp_end,3=HOM_count,4=HET_count,5=oddHET_count]
					if C_het==False:
						fragments[current_fragment-1][3] += 1
					else:
						if C_oddHet==False:
							fragments[current_fragment-1][4] += 1
						else:
							fragments[current_fragment-1][5] += 1

				# If current coordinate is at end of a fragment, update fragment_found to 0.
				if pos_point == fragments[current_fragment-1][2]:
					fragment_found = 0

		# Reset old_chr/last_fragment to current coordinate chromosome before moving to next line in pileup. 
		old_chr       = chr
		last_fragment = current_fragment

print "### ", time.clock() - t1, "seconds to parse project SNP data."
t2 = time.clock()
print '### Number of fragments = ' + str(numFragments)
print '### Data from each fragment: [chrNum, bpStart, bpEnd, Max, Ave, Length]'


#============================================================================================================
# Code section to output information about read average per standard bin genome fragment.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\n\tOutputting LOH status counts of standard-bin fragmented genome.\n")
print '### chr_num\tbp_start\tbp_end\tHOM_count\tHET_count\toddHET_count'
for fragment in range(1,numFragments):
	# Output a line for each fragment.
	#     fragments[fragment-1] = [chr_num,bp_start,bp_end, GC_ratio]
	#     0) chr_num
	#     1) bp_start
	#     2) bp_end
	#     3) HOM_count
	#     4) HET_count
	#     5) oddHET_count
	chr_num      = fragments[fragment-1][0]
	bp_start     = fragments[fragment-1][1]
	bp_end       = fragments[fragment-1][2]
	HOM_count    = fragments[fragment-1][3]
	HET_count    = fragments[fragment-1][4]
	oddHET_count = fragments[fragment-1][5]

	print str(chr_num) + '\t' + str(bp_start) + '\t' + str(bp_end) + '\t' + str(HOM_count) + '\t' + str(HET_count) + '\t' + str(oddHET_count)

#------------------------------------------------------------------------------------------------------------
# End of code section to output information about fragments. 
#============================================================================================================

print "### ", time.clock() - t1, "seconds to output basic stats of each restriction fragment."
print "### ", time.clock() - t0, "seconds to complete processing of fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\n\tTime to process = " + str(time.clock()-t0) )
	myfile.write("\n\t* 'py/dataset_process_for_CNV_analysis.py' completed. *\n")
