def processLine(entry_line):
	global chrNums
	global chrName
	global chrCount

	# Process 'SNP_CNV_v1.txt' file line.
	# example lines from SNP_CNV pileup file:
	#       chromosome               coord   total   ref   A    T    G    C
	#       ChrA_C_glabrata_CBS138   45      36      C     0    0    0    36
	#       ChrA_C_glabrata_CBS138   46      37      T     0    37   0    0
	#       ChrA_C_glabrata_CBS138   47      38      A     38   0    0    0
	#       ChrA_C_glabrata_CBS138   48      39      A     39   0    0    0

	parent_line = string.strip(entry_line)
	parent_line = parent_line.split('\t')
	P_chr_name  = parent_line[0]   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	P_position  = parent_line[1]   # chr position of bp.      : 2286371
	P_countTot  = parent_line[2]   # total count at bp.       : 101
	P_refBase   = parent_line[3]   # reference base at bp.    : T
	P_countA    = parent_line[4]   # count of A.              : 100
	P_countT    = parent_line[5]   # count of T.              : 0
	P_countG    = parent_line[6]   # count of G.              : 0
	P_countC    = parent_line[7]   # count of C.              : 1
	# Determine chrID associated with chromosome name.
	P_chr = 0
	for x in range(0,chrCount):
		if (chrNums[x] != 0):
			if chrName[x] == P_chr_name:
				P_chr = x+1
	P_chrName = chrName[P_chr-1]

	return P_chr,P_chrName,P_position,P_countA,P_countT,P_countG,P_countC


###
### Preprocesses parental 'SNP_CNV_v1.txt' files that has been copied into the hapmap folder
###     as "SNPdata_parent1.txt" and "SNPdata_parent2.txt"
### Final format only has data lines for loci which are homozygous in each parent, meaning an allelic
###     ratio of greater than 0.5 for major allele, and different between each parent.
### Output collumns:
###     0) chr_name      : Chromosome name string.
###     1) bp_coordinate : bp coordinate along chromosome.
###     2) baseCall_1    : alphabetical first call of heterozgous locus.
###     3) baseCall_2    : alphabetical second call of heterozygous locus.
### Uses genome definition files to only output data lines for chromosomes of interest.
###
### Inputfile1 & Inputfile2 are 'SNP_CNV_v1.txt' files for the two parent strains.
###

import string, sys, time

genome       = sys.argv[ 1]
genomeUser   = sys.argv[ 2]
project1     = sys.argv[ 3]
project1User = sys.argv[ 4]
project2     = sys.argv[ 5]
project2User = sys.argv[ 6]
hapmap       = sys.argv[ 7]
hapmapUser   = sys.argv[ 8]
main_dir     = sys.argv[ 9]
runMode      = sys.argv[10]

logName     = main_dir+"users/"+hapmapUser+"/hapmaps/"+hapmap+"/process_log.txt"
inputFile1  = main_dir+"users/"+hapmapUser+"/hapmaps/"+hapmap+"/SNPdata_parent1.txt"
inputFile2  = main_dir+"users/"+hapmapUser+"/hapmaps/"+hapmap+"/SNPdata_parent2.txt"
with open(logName, "a") as myfile:
	myfile.write("|\trunMode = 'hapmap'\n")
	myfile.write("|\t    Building hapmap from '"+project1+"' and '"+project2+"',\n")
	myfile.write("|\t        defining homolog 'a' and 'b', respespectively.\n")

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("*==================================================*\n")
	myfile.write("| Log of 'py/hapmap.preprocess_parent.py'          |\n")
	myfile.write("*--------------------------------------------------*\n")


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

#============================================================================================================
# Process 'preprocessed_SNPs.txt' file for projectParent to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tProcessing parent 'preprocessed_SNPs.txt' file -> het loci.\n")

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
	myfile.write("|\tDetermining number of chromosomes of interest in genome.\n")

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
				myfile.write("|\t\t" + str(chr_num) + " : " + chr_name + " = " + chr_nameShort + "\n")
figureDefinitionFile.close()

# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount
with open(logName, "a") as myfile:
	myfile.write("|\t\tMax chr string : "+str(chrCount)+"\n")
#............................................................................................................

with open(logName, "a") as myfile:
	myfile.write("|\tOpen parent 'putative_SNPs_v4.txt' file.\n")

#............................................................................................................

count            = 0
old_chr          = 0
fragment_found   = 0
last_fragment    = 0
current_fragment = 0
log_count        = 0
log_offset       = 0

print '### Chromosomes of interest : '
for x in range(0,chrCount):
	if (chrNums[x] != 0):
		print '### \t' + str(x+1) + ' : ' + str(chrName[x])

with open(logName, "a") as myfile:
	myfile.write("|\tGathering read coverage data for each fragment.\n")

# Open dataset 'putative_CNVs_v1.txt' file.
data1 = open(inputFile1,"r")
data2 = open(inputFile2,"r")

print '### Data lines for each het locus in parent : [chromosome_name, bp_coordinate, allele_1, allele_2]'

# Process 'SNP_CNV_v1.txt' file for both parents, line by line... while checking for missing data.
line1 = data1.readline()
line2 = data2.readline()
error_endOfFile = False
while (error_endOfFile == False):
	#	parent1_line = string.strip(line1)
	#	parent1_line = parent1_line.split('\t')
	#	P1_chr_name  = parent1_line[0]   # chr name of bp.			: Ca21chrR_C_albicans_SC5314
	#	P1_position  = parent1_line[1]   # chr position of bp.		: 2286371
	#	P1_countTot  = parent1_line[2]   # total count at bp.       : 101
	#	P1_refBase   = parent1_line[3]   # reference base at bp.	: T
	#	P1_countA    = parent1_line[4]   # count of A.				: 100
	#	P1_countT    = parent1_line[5]   # count of T.				: 0
	#	P1_countG    = parent1_line[6]   # count of G.				: 0
	#	P1_countC    = parent1_line[7]   # count of C.				: 1
	#	# Determine chrID associated with chromosome name.
	#	P1_chrID = 0
	#	for x in range(0,chrCount):
	#		if (chrNums[x] != 0):
	#			if chrName[x] == P1_chr_name:
	#				P1_chr = x+1

	P1_chrID,P1_chrName,P1_position,P1_countA,P1_countT,P1_countG,P1_countC = processLine(line1)
	P2_chrID,P2_chrName,P2_position,P2_countA,P2_countT,P2_countG,P2_countC = processLine(line2)

#	print 'A ' + P1_chrName +'\t'+ str(P1_position) +'\t'+ P1_countA +':'+ P1_countT +':'+ P1_countG +':'+ P1_countC
#	print 'A ' + P2_chrName +'\t'+ str(P2_position) +'\t'+ P2_countA +':'+ P2_countT +':'+ P2_countG +':'+ P2_countC

	## Ensure that we're at the same chromosome in each file.
	#  Files are sorted by chromosome, with increasing chromosome IDs.
	#  If file load line hits the end of the file, then break out of the loops, as we're
	#      done with file processing.
#	print 'place 1'
	while (P1_chrID <> P2_chrID):
#		print 'place 2'
		while (P1_chrID < P2_chrID):
	#		print 'place 3'
			# load a line from inputFile1.
			line1 = data1.readline()
			if not line1: # EOF 1
				error_endOfFile = True
				break
			P1_chrID,P1_chrName,P1_position,P1_countA,P1_countT,P1_countG,P1_countC = processLine(line1)
		while (P2_chrID < P1_chrID):
	#		print 'place 4'
			# load a line from inputFile2.
			line2 = data2.readline()
			if not line2: # EOF 2
				loop_eror = True
				break
			P2_chrID,P2_chrName,P2_position,P2_countA,P2_countT,P2_countG,P2_countC = processLine(line2)
		if error_endOfFile:
			break
	# End of while loop that registers file pointer to chromosome.
	if error_endOfFile:
		break
	## Ensure we're at the same bp coordinate in each file.
	#  Files are sorted by coordinate, with increaseing coordinate bps.
	#  If file load line hits the end of the file, then break out of the loops, as we're
	#      done with file processing.
	#  If file load line hits a new chromosome, then break
	error_endOfFile = False
	error_endOfChr  = False
	# print 'place 5'
	while (P1_position <> P2_position):
	#	print 'place 6'
		while (P1_position < P2_position):
	#		print 'place 7'
			# load a line from inputFile1.
			line1 = data1.readline()
			if not line1: # endOfFile 1
				error_endOfFile = True
				break;
			P1_chrID,P1_chrName,P1_position,P1_countA,P1_countT,P1_countG,P1_countC = processLine(line1)
			if (P1_chrID <> P2_chrID): # endOfChr 1
				error_endOfChr = True
				break;
		while (P2_position < P1_position):
	#		print 'place 8'
			# load a line from inputFile2.
			line2 = data2.readline()
			if not line2: # endOfFile 2
				loop_eror = True
				break
			P2_chrID,P2_chrName,P2_position,P2_countA,P2_countT,P2_countG,P2_countC = processLine(line2)
			if (P1_chrID <> P2_chrID): # endOfChr 2
				error_endOfChr = True
				break
		if error_endOfFile:
			break
		if error_endOfChr:
			break
	# End of while loop that registers file pointer to coordinate along chromosome.
	# print 'place 9'
	if error_endOfFile:
		break

	# At this time, line and line2 shuold point to the same chromosome and coordinate, or the
	# file processing loop will have exited because the locus was not found in one or the other
	# input data file.
	# Between the two data files, this locus may be het or hom...  but we can't tell, so it is
	# not added to the output heterozygous loci list.

	P1_list = [int(float(P1_countA)), int(float(P1_countT)), int(float(P1_countG)), int(float(P1_countC))]
	if (sum(P1_list) == 0):
		P1_allelicRatio = 0
	else:
		P1_allelicRatio = max(P1_list)/float(sum(P1_list))
	P2_list = [int(float(P2_countA)), int(float(P2_countT)), int(float(P2_countG)), int(float(P2_countC))]
	if (sum(P2_list) == 0):
		P2_allelicRatio = 0
	else:
		P2_allelicRatio = max(P2_list)/float(sum(P2_list))
	if not ((P1_allelicRatio < 0.5) or (P2_allelicRatio < 0.5)):
	#	print 'place 10'
		# At this time, any loci which appear to be heterozygous in either haploid parent will have
		# been excluded. Such sites are problematic and likely due to duplicated genes/regions within
		# the genome. Such 'SNPs' aren't useful in defining the haplotype map.

		P1_alleleData    = [(int(float(P1_countA)),'A'), (int(float(P1_countT)),'T'), (int(float(P1_countG)),'G'), (int(float(P1_countC)),'C')]
		P1_sortedAlleles = sorted(P1_alleleData, key=lambda alleleDatum: alleleDatum[0]) # sort alleles by copy number.
		P1_majorAllele   = P1_sortedAlleles[3][1]
		P2_alleleData    = [(int(float(P2_countA)),'A'), (int(float(P2_countT)),'T'), (int(float(P2_countG)),'G'), (int(float(P2_countC)),'C')]
		P2_sortedAlleles = sorted(P2_alleleData, key=lambda alleleDatum: alleleDatum[0]) # sort alleles by copy number.
		P2_majorAllele   = P2_sortedAlleles[3][1]
		if not (P1_majorAllele == P2_majorAllele):
	#		print 'place 11'
			# At this point, any loci which are homozygous and identical in both haploid parents will have
			# been excluded. Such sites are of no use in defining the haplotype map differences between the
			# haploid parents.

			# If the locus remains at this point, we add a line for it to the output pileup file.
			# Output new pileup line
			#	0) chr_name      : Chromosome name string.
			#	1) bp_coordinate : bp coordinate along chromosome.
			#	2) baseCall_1    : allele in haplotype 'a', from dataset1.
			#	3) baseCall_2    : allele in haplotype 'b', from dataset2.
			if P1_chrID!=0:
				print P1_chrName + '\t' + str(P1_position) + '\t' + P1_majorAllele + '\t' + P2_majorAllele + '\t0'
				# print P1_chrName +'\t'+ str(P1_position) +'\t['+ str(P1_countA) +','+ str(P1_countT) +','+ str(P1_countG) +','+ str(P1_countC) +']:['+ str(P2_countA) +','+ str(P2_countT) +','+ str(P2_countG) +','+ str(P2_countC) +']\t'+ P1_majorAllele +'\t'+ P2_majorAllele
				# print '{'+ P1_sortedAlleles[0][1] +','+ P1_sortedAlleles[1][1] +','+ P1_sortedAlleles[2][1] +','+ P1_sortedAlleles[3][1] +'}:{'+ P2_sortedAlleles[0][1] +','+ P2_sortedAlleles[1][1] +','+ P2_sortedAlleles[2][1] +','+ P2_sortedAlleles[3][1] +'}'
				# Final column of '0' is added to indicate no need to correct phasing.
	line1 = data1.readline()
	if not line1: # EOF 1
		error_endOfFile = True
		break
	line2 = data2.readline()
	if not line2: # EOF 2
		error_endOfFile = True
		break

#------------------------------------------------------------------------------------------------------------
# End of code section to output information about parental heterozygous loci. 
#============================================================================================================

print '### End of preprocessed parental SNP data.'

with open(logName, "a") as myfile:
	myfile.write("|\tTime to process = " + str(time.clock()-t0) + "\n")
	myfile.write("*--------------------------------------------------*\n")
	myfile.write("| End of 'py/hapmap.preprocess_parent.py'          |\n")
	myfile.write("*=================================================-*\n")
