# Input Arguments
#   1) user     : 'darren'
#   2) project  : 'test_Ca'
#   3) genome   : 'Candida_albicans_SC5314_etc'
#   4) main_dir : '/home/bermanj/shared/links/'
#	5) logName  :
#
# Process input files:
#	1) Raw SNP data                 : $workingDir"users/"$user"/projects/"$project"/SNP_CNV_v1.txt".
#	2) FASTA file name              : $workingDir"users/default/genomes/default/reference.txt",
#	                               or $workingDir"users/"$user"/genomes/default/reference.txt" as $FastaName.
#	3) Coordinates of standard bins : $workingDir"users/default/genomes/"$genome"/"$FastaName"."$restrictionEnzymes".fasta",
#	                               or $workingDir"users/"$user"/genomes/"$genome"/"$FastaName"."$restrictionEnzymes".fasta".

# Generate output file:
#	1) a simplified pileup file containing data about parental het loci in child dataset per RADseq fragment.
#		0) chr_num  : Numerical chromosome identifier, defined for each genome in "figure_details.txt".
#		1) bp_start : Start bp coordinate along chromosome.
#		2) bp_end   : End bp coordinate along chromosome.
#		3) list of phased data ratios.
#		4) list of unphased data ratios.
#		5) list of phased data coordinates.
#		6) list of unphased data coordinates.
#		Comment lines in output begin with '#'.
#

import string, sys, re, time, os
genome             = sys.argv[ 1]
genomeUser         = sys.argv[ 2]
hapmap             = sys.argv[ 3]	# 'hapmap' or 'parent project'...            both will contain a 'putative_SNPs_v4.txt'/'SNPdata_parent.txt' format file.
hapmapUser         = sys.argv[ 4]	# 'hapmap user' or 'parent project user'... 
project            = sys.argv[ 5]	# 'child project'.
projectUser        = sys.argv[ 6]	# 'child project user'.
main_dir           = sys.argv[ 7]
restrictionEnzymes = sys.argv[ 8]
logName            = sys.argv[ 9]
runMode            = sys.argv[10]  # 'hapmap' or 'LOH' modes.

with open(logName, "a") as myfile:
	myfile.write("*---------------------------------------------------------*\n")
	myfile.write("| Log of 'dataset_process_for_SNP_analysis.ddRADseq.py'   |\n")
	myfile.write("*---------------------------------------------------------*\n")

# Figure out if input name 'hapmap' corresponds to a project or an actual hapmap.
if (runMode == 'hapmap'):
	parentDatafile  = main_dir+"users/"+hapmapUser+"/hapmaps/"+hapmap+"/SNPdata_parent.txt"
	with open(logName, "a") as myfile:
		myfile.write("|\trunMode = 'hapmap'\n")
		myfile.write("|\t    Comparing project '"+project+"' to hapmap '"+hapmap+"'.\n")
elif (runMode == 'LOH'):
	parentDatafile = main_dir+"users/"+projectUser+"/projects/"+project+"/SNPdata_parent.txt"
	with open(logName, "a") as myfile:
		myfile.write("|\trunMode = 'LOH'\n")
		myfile.write("|\t    Comparing project '"+project+"' to parent project '"+hapmap+"'.\n")

childDatafile = main_dir+"users/"+projectUser+"/projects/"+project+"/SNP_CNV_v1.txt"

t0 = time.clock()


#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+genomeUser+"/genomes/"+genome+"/"

with open(logName, "a") as myfile:
	myfile.write("|\tProcessing standard bin fragmented genome file.\n")


#============================================================================================================
# Load FastaName from 'reference.txt' for genome in use.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tIdentifying name of reference FASTA file.\n")
reference_file = genomeDirectory + '/reference.txt'
refFile        = open(reference_file,'r')
FastaName      = refFile.read().strip()
refFile.close()
FastaName      = FastaName.replace(".fasta", "")
with open(logName, "a") as myfile:
	myfile.write("|\t\t'"+FastaName+"'\n")


#============================================================================================================
# Process restriction-digested genome file.
#------------------------------------------------------------------------------------------------------------
# Example FASTQ header line.
#     >Ca_a.chr1 (9638..10115) (478bp) [*]
with open(logName, "a") as myfile:
    myfile.write("\t\t\t|\tLoading ddRADseq fragmented genome FASTA file.\n")
# Open restriction-digested genome FASTA file.
ddRADseq_FASTA_file = genomeDirectory + FastaName + "." + restrictionEnzymes + ".fasta"
ddRADseq_FASTA_data = open(ddRADseq_FASTA_file,'r')
# Setup array and counter for tracking fragment definition data.
fragments        = []
fragment_counter = 0
## Process digested FASTA genome file, line by line.
while True:
	# Line pairs have the following structure.
	#    >Ca_a.chr1 (9638..10115) (478bp) [*]
	#    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
	line1 = ddRADseq_FASTA_data.readline()
	line2 = ddRADseq_FASTA_data.readline()
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
				coordinates        = bp_coordinate_string.replace('(','').replace(')','').replace('..',' ').split()
				bp_start           = int(float(coordinates[0]))
				bp_end             = int(float(coordinates[1]))

				phasedData         = '(' # start of string for phased data.
				unphasedData       = '(' # start of string for unphased data.
				phasedCoordinate   = '(' # start of string for phased data coordinates.
				unphasedCoordinate = '(' # start of string for unphased data coordinaes.
				fragments.append([chr_num,bp_start,bp_end, phasedData,unphasedData,phasedCoordinate,unphasedCoordinate])
				#print '###\tfragment[' + str(fragment_counter) + '] = [' + str(chr_num) + ', ' + str(bp_start) + ', ' + str(bp_end) + ', ' + ']'
				fragment_counter += 1
ddRADseq_FASTA_data.close()
with open(logName, "a") as myfile:
	myfile.write("|\tThe ddRADseq fragmented genome FASTA file has been loaded.\n")
# Put fragment counter into a general use variable.
numFragments = fragment_counter
#------------------------------------------------------------------------------------------------------------
# End of code section to parse standard-bin fragments file.
#============================================================================================================


print "### ", time.clock() - t0, "seconds to parse restriction fragments from digested genome."
t1 = time.clock()
print "### Starting read count data processing."


#============================================================================================================
# Process hapmap 'SNPdata_parent.txt' file to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tProcessing 'SNPdata_parent.txt' file -> het loci.\n")
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
# Determine the number of chromosomes of interest in genome.
with open(logName, "a") as myfile:
	myfile.write("|\tDetermining number of chromosomes of interest in genome.\n")
chrName_maxcount = 0
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line))
	chr_num = line_parts[0]
	if chr_num.isdigit():
		chr_num    = int(float(line_parts[0]))
		chr_use    = int(float(line_parts[1]))
		chr_label  = line_parts[2]
		chr_name   = line_parts[3]
		if chr_num > chrName_maxcount:
			chrName_maxcount = chr_num
figureDefinitionFile.close()
# Pre-allocate chrName_array
chrName = []
for x in range(0, chrName_maxcount+10):
	chrName.append([])
with open(logName, "a") as myfile:
	myfile.write("|\tGathering name strings for chromosomes.\n")
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
		chr_use                        = int(float(line_parts[1]))
		chr_label                      = line_parts[2]
		chrLabels.append(chr_label);
		chr_name                       = line_parts[3]
		chrNames.append(chr_name);
		chr_nameShort                  = chr_label
		chrShorts.append(chr_nameShort);
		chrName[chrCounter] = chr_name
		with open(logName, "a") as myfile:
			myfile.write("|\t\t" + str(chr_num) + " : " + chr_name + " = " + chr_nameShort + "\n")
		chrCounter += 1
figureDefinitionFile.close()
# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount
with open(logName, "a") as myfile:
	myfile.write("|\t\tMax chr string : "+str(chrCount)+"\n")
	myfile.write("|\tProcessing datasetParent 'SNPdata_parent.txt' file => het loci.\n")
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
	myfile.write("|\tGathering read coverage data for each fragment.\n")
	myfile.write("|\tparent data file = '"+parentDatafile+"'\n")
	myfile.write("|\tchild data file  = '"+childDatafile+"'\n")
# Open hapmap 'SNPdata_parent.txt' file, which only contains lines for heterozygous SNP loci in the parent dataset.
print '### parentDatafile = ' + parentDatafile
data         = open(parentDatafile,'r')
searchTarget = open(childDatafile,'r')
childLine    = ''

# Process hapmap 'SNPdata_parent.txt' file, line by line.
for line in data:
	#===============================================================================================================
	# Process the line from the parent 'SNPdata_parent.txt' line, describing one heterozgyous locus.
	#--------------------------------------------------------------------------------------------------------------- 
	# example lines from file:
	#       chromosome                   coord   allele1   allele2   hapmapEntry1
	#       Ca21chr1_C_albicans_SC5314   3706    T         C         1
	if line[0] != "#":
		count += 1
		parentLine    = string.strip(line)
		parentLine    = parentLine.split('\t')
		P_chr_name    = parentLine[0]        # chr name of bp.		: Ca21chrR_C_albicans_SC5314
		P_position    = int(parentLine[1])   # chr position of bp.	: 2286371
		P_allele1     = parentLine[2]        # allele 1.			: T
		P_allele2     = parentLine[3]        # allele 2.			: A
		phasingData_1 = parentLine[4:]       # list of phasing data points.   remove any non-useful data points (10,11,12).   resulting in, 0s and 1s...  sum phase info = round(sum(list)/len(list))
		phasingData_2 = [x for x in phasingData_1 if int(x) != 10]   # remove no phase condition '10' = het coordinate not found in dataset.
		phasingData_3 = [x for x in phasingData_2 if int(x) != 11]   # remove no phase condition '11' = het coordinate not associated with LOH fragment definition.
		phasingData   = [x for x in phasingData_3 if int(x) != 12]   # remove no phase condition '12' = het coordinate allele not in hapmap.
		# Determine summary phase call for locus, used to apply counts to proper output columns.
		if len(phasingData) == 0:
			# no phasing data available.
			phaseCall = 10 # homolog undefined.
		else:
			# at least one phasing data point.
			phasingData = [int(i) for i in phasingData]
			dataLength  = len(phasingData)
			dataSum     = sum(phasingData)
			dataAve     = round(dataSum/float(dataLength))
			phaseCall   = int(dataAve) #  [0,1] for homologs 'a' and 'b'.
		# Identify which chromosome this data point corresponds to.
		P_chr = 0
		for x in range(0,chrCount):
			if (chrNums[x] != 0):
				if chrName[x] == P_chr_name:
					P_chr = x+1
		####################################################
		# All hapmap dataset coordinates are heterozygous. #
		####################################################
		# Find standard bin containing parent SNP locus.
		if P_chr != 0:
			# Reset for each new chromosome examined.
			if old_chr != P_chr:
				if log_offset != 0:
					log_offset_string = " "*((log_offset)%100)
					with open(logName, "a") as myfile:
						myfile.write("\n|\t\t" + log_offset_string)
				count            = 1
				fragment_found   = 0
				current_fragment = 0
			# If (fragment_found == 0), look to see if current coordinate matches any defined fragments.
			if fragment_found == 0:
				for frag in range(current_fragment,numFragments):
					# Check if current coordinate is consistent with this fragment : fragments[frag-1] = [chr_num,bp_start,bp_end, data_count,data_max,ave_read_count, repetDataCount,repetMax,repetAve]
					if P_chr == fragments[frag-1][0] and P_position >= fragments[frag-1][1] and P_position <= fragments[frag-1][2]:
						fragment_found   = 1
						current_fragment = frag
						break
			#print str(numFragments)+":"+str(current_fragment)
			# If (fragment_found == 1), add current bp coordinate data to fragment data.
			if fragment_found == 1:
				# display status updates to log file.
				if (count%50000 == 0):
					log_count  += 1
					log_offset += 1
					if (log_count%100) == 0:
						if (log_count-1) == 0:
							with open(logName, "a") as myfile:
								myfile.write("\n|\t")
						else:
							with open(logName, "a") as myfile:
								myfile.write(" " + str(log_count-1))
					with open(logName, "a") as myfile:
						myfile.write(".")
				# If current coordinate is at/after end of a fragment, update fragment_found to 0.
				if P_position >= fragments[current_fragment-1][2]:
					fragment_found = 0
			#===============================================================================================================
			# Find line in dataset corresponding to SNP from hapmap dataset 'SNPdata_parent.txt' file.
			#--------------------------------------------------------------------------------------------------------------- 
			searchString = P_chr_name+"\t"+str(P_position)
			result = ""
			## Barely functional, very slow slow equivalent of GREP, returning only one result line.
			##    Very slow because it restarts search position for each search.
			##    This can be improved dramatically because the search target file is already sorted.
			#for lineChild in open(inputFile2,'r'):
			#   if searchString in lineChild:
			#       result = lineChild
			#       break
			# Fast line search, relies upon presorted search target file.
			if len(childLine) == 0:
				childLine       = searchTarget.readline()
			childLine_parts     = childLine.split('\t')
			C_chr_name          = childLine_parts[0]
			C_position          = int(childLine_parts[1])
			C_chr = 0
			for x in range(0,chrCount):
				if (chrNums[x] != 0):
					if chrName[x] == C_chr_name:
						C_chr = x+1
			#print "1|P:C "+str(P_chr)+":"+str(C_chr)+" "+str(P_position)+":"+str(C_position)+"|"
			while P_chr > C_chr:    # WORKING: this section jumps through the child lines of chromosomes with no parent lines.
				childLine       = searchTarget.readline()
				if len(childLine) > 0:
					childLine       = string.strip(childLine)
					childLine_parts = childLine.split('\t')
					C_chr_name      = childLine_parts[0]
					C_position      = int(childLine_parts[1])
					for x in range(0,chrCount):
						if (chrNums[x] != 0):
							if chrName[x] == C_chr_name:
								C_chr = x+1
				else:
					C_chr = P_chr
			#print "2|P:C "+str(P_chr)+":"+str(C_chr)+" "+str(P_position)+":"+str(C_position)+"|"
			while P_chr == C_chr and P_position > C_position:   # WORKING: this section jumps through the child lines until the correct chromosome and coordinate is reached.
				childLine       = searchTarget.readline()
				if len(childLine) > 0:
					childLine       = string.strip(childLine)
					childLine_parts = childLine.split('\t')
					C_chr_name      = childLine_parts[0]
					C_position      = int(childLine_parts[1])
					for x in range(0,chrCount):
						if (chrNums[x] != 0):
							if chrName[x] == C_chr_name:
								C_chr = x+1
				else:
					C_position = P_position
			#print "3|P:C "+str(P_chr)+":"+str(C_chr)+" "+str(P_position)+":"+str(C_position)+"|"
			if P_chr > C_chr:
				result          = ""
			elif P_chr == C_chr and P_position > C_position:
				result          = ""
			elif P_chr == C_chr and P_position < C_position:
				result          = ""
			elif P_chr == C_chr and P_position == C_position:
				result          = childLine
			#print "4|"+childLine+"|"
			C_ratio        = 0.0
			if result == "":
				# locus is not found in dataset, no contribution to SNP interpretations.
				# print '# locus not found in dataset.'
				C_ratio        = 0.0
				C_valid        = 0
			else:
				C_chr_name     = childLine_parts[0]   # chr_name_string
				C_position     = childLine_parts[1]   # bp_coordinate
				C_countTot     = childLine_parts[2]   # total_reads
				C_refBase      = childLine_parts[3]   # reference_base
				C_countA       = childLine_parts[4]   # A reads
				C_countT       = childLine_parts[5]   # T reads
				C_countG       = childLine_parts[6]   # G reads
				C_countC       = childLine_parts[7]   # C reads
				C_counts       = [int(float(C_countA)), int(float(C_countT)), int(float(C_countG)), int(float(C_countC))]
				# Child locus should be homozygous, or at least unbalanced, so figure out which allele is most common.
				alleleData     = [(int(float(C_countA)), 'A'), (int(float(C_countT)), 'T'), (int(float(C_countG)), 'G'), (int(float(C_countC)), 'C')]
				sortedAlleles  = sorted(alleleData, key=lambda alleleDatum: alleleDatum[0]) # sort alleles by copy number.

				C_count1       = [item for item in alleleData if item[1] == P_allele1][0][0]
				C_count2       = [item for item in alleleData if item[1] == P_allele2][0][0]
				C_ratio = 1
				if sum(C_counts) == 0:
					C_ratio    = 0.0
					C_valid    = 0
				else:
					if phaseCall == 0:
						C_ratio    = float(C_count1)/float(sum(C_counts))
					elif phaseCall == 1:
						C_ratio    = float(C_count2)/float(sum(C_counts))
					else:
						C_ratio    = float(C_count1)/float(sum(C_counts))
					C_valid    = 1
				#print '# locus found : '+P_allele1+'/'+P_allele2+' : '+str(C_counts)+' : '+str(C_count1)+'/'+str(C_count2)+'/'+str(sum(C_counts))+' : '+str(C_ratio)

			C_chr = 0
			for x in range(0,chrCount):
				if (chrNums[x] != 0):
					if chrName[x] == P_chr_name:
						C_chr = x+1

			#===============================================================================================================
			# Add allelic ratio data to standard-bin fragment data
			#---------------------------------------------------------------------------------------------------------------
			if fragment_found == 1:
				if  C_valid == 1:
					# standard-bin fragments structure:   'current_fragment-1'
					#   fragments    = [chr_num,bp_start,bp_end,phasedData,unphasedData,phasedCoordinate,unphasedCoordinate]
					if phaseCall == 0 or phaseCall == 1:
						# correctly phased data, or incorrectly phased data.
						#    the corrected C_ratio is calculated in the above section.
						fragments[current_fragment-1][3] = fragments[current_fragment-1][3]+str(C_ratio)+','
						fragments[current_fragment-1][5] = fragments[current_fragment-1][5]+str(C_position)+','
					else:
						# unphased data.
						fragments[current_fragment-1][4] = fragments[current_fragment-1][4]+str(C_ratio)+','
						fragments[current_fragment-1][6] = fragments[current_fragment-1][6]+str(C_position)+','
			#print "|"+str(current_fragment)+" "+str(C_ratio)

			#===============================================================================================================
			# Output log file status updates.
			#---------------------------------------------------------------------------------------------------------------
			if old_chr != P_chr:
				print '### chr change : ' + str(old_chr) + ' -> ' + str(P_chr)
				with open(logName, "a") as myfile:
					myfile.write("\n|\t    " + str(old_chr) + " -> " + str(P_chr) + " = " + P_chr_name + "\n")
					myfile.write(  "|\t1........01........01........01........01........01........01........01........01........01........0")
			# Reset for each new chromosome examined.
			if old_chr != P_chr:
				if log_offset != 0:
					log_offset_string = " "*(log_offset%100)
					with open(logName, "a") as myfile:
						myfile.write("\n|\t" + log_offset_string)
			# If (fragment_found == 1), add current bp coordinate data to fragment data.
			if fragment_found == 1:
				# display status updates to log file.
				if last_fragment != current_fragment:
					log_count  += 1
					log_offset += 1
					if ((log_count-1)%100) == 0:
						if (log_count-1) == 0:
							with open(logName, "a") as myfile:
								myfile.write("\n|\t")
						else:
							with open(logName, "a") as myfile:
								myfile.write(" " + str(log_count-1) + "\n|\t")
					with open(logName, "a") as myfile:
						myfile.write(".")

		#===============================================================================================================
		# Update tracking variables to current coordinates before moving to next line in 'SNPdata_parent.txt' file.
		#---------------------------------------------------------------------------------------------------------------
		old_chr       = P_chr
		last_fragment = current_fragment
		# If current coordinate is at end of a fragment, update fragment_found to 0.
		if P_position == fragments[current_fragment-1][2]:
			fragment_found = 0

print "### ", time.clock() - t1, "seconds to parse project SNP data."
t2 = time.clock()
print '### Number of fragments = ' + str(numFragments)
print '### Data from each fragment: [chrNum, bpStart, bpEnd, Max, Ave, Length]'


#============================================================================================================
# Code section to output data per ddRADseq fragment.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\n|\tOutputting LOH status counts of standard-bin fragmented genome.\n")
print '### chr_num\tbp_start\tbp_end\tHOM_count\tHET_count\toddHET_count'
for fragment in range(1,numFragments):
	# Output a line for each fragment.
	#     fragments[fragment-1] = [chr_num,bp_start,bp_end,phasedData,unphasedData]
	#     0) chr_num
	#     1) bp_start
	#     2) bp_end
	#     3) list of phased data ratios.
	#     4) list of unphased data ratios.
	#     5) list of phased data coordinates.
	#     6) list of unphased data coordinates.
	chrNum_string           =  str(fragments[fragment-1][0])
	bpStart_string          =  str(fragments[fragment-1][1])
	bpEnd_string            =  str(fragments[fragment-1][2])
	phasedData_string       =  fragments[fragment-1][3]
	if len(phasedData_string) > 1:
		phasedData_string   =  phasedData_string[:-1]
	phasedData_string       += ")"
	unphasedData_string     =  fragments[fragment-1][4]
	if len(unphasedData_string) > 1:
		unphasedData_string       =  unphasedData_string[:-1]
	unphasedData_string           += ")"
	phasedCoordinate_string       =  fragments[fragment-1][5]
	if len(phasedCoordinate_string) > 1:
		phasedCoordinate_string   =  phasedCoordinate_string[:-1]
	phasedCoordinate_string       += ")"
	unphasedCoordinate_string     =  fragments[fragment-1][6]
	if len(unphasedCoordinate_string) > 1:
		unphasedCoordinate_string =  unphasedCoordinate_string[:-1]
	unphasedCoordinate_string     += ")"

	print chrNum_string+'\t'+bpStart_string+'\t'+bpEnd_string+'\t'+phasedData_string+'\t'+unphasedData_string+'\t'+phasedCoordinate_string+'\t'+unphasedCoordinate_string

#------------------------------------------------------------------------------------------------------------
# End of code section to output information about fragments. 
#============================================================================================================

print "### ", time.clock() - t1, "seconds to output basic stats of each restriction fragment."
print "### ", time.clock() - t0, "seconds to complete processing of fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("|\tTime to process = " + str(time.clock()-t0) +"\n")
	myfile.write("*--------------------------------------------------------------*\n")
	myfile.write("| 'py/dataset_process_for_SNP_analysis.ddRADseq.py' completed. |\n")
	myfile.write("*--------------------------------------------------------------*\n")
