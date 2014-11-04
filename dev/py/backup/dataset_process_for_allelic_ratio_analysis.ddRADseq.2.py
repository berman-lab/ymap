def processLine(entry_line):
	global chrNums
	global chrName
	global chrCount

	# example lines from 'putative_SNPs_v4.txt' pileup file:
	#       chromosome               coord   ref   A    T    G    C
	#       ChrA_C_glabrata_CBS138   45      C     0    0    0    36
	#       ChrA_C_glabrata_CBS138   46      T     0    37   0    0
	#       ChrA_C_glabrata_CBS138   47      A     38   0    0    0
	#       ChrA_C_glabrata_CBS138   48      A     39   0    0    0

	parent_line = string.strip(entry_line)
	parent_line = parent_line.split('\t')
	P_chr_name  = parent_line[0]   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	P_position  = parent_line[1]   # chr position of bp.      : 2286371
	P_refBase   = parent_line[2]   # reference base at bp.    : T
	P_countA    = parent_line[3]   # count of A.              : 100
	P_countT    = parent_line[4]   # count of T.              : 0
	P_countG    = parent_line[5]   # count of G.              : 0
	P_countC    = parent_line[6]   # count of C.              : 1
	# Determine chrID associated with chromosome name.
	P_chr = 0
	for x in range(0,chrCount):
		if (chrNums[x] != 0):
			if chrName[x] == P_chr_name:
				P_chr = x+1
	P_chrName = chrName[P_chr-1]

	return P_chr,P_chrName,P_position,P_countA,P_countT,P_countG,P_countC

def processLine_CNV_SNP(entry_line):
	global chrNums
	global chrName
	global chrCount

	parent_line = string.strip(entry_line)
	parent_line = parent_line.split('\t')
	P_chr_name  = parent_line[0]   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	P_position  = parent_line[1]   # chr position of bp.      : 2286371
	P_countTot  = parent_line[2]
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
### Output collumns:
###     0) chr_name        : Chromosome name string.
###     1) bp_coordinate   : bp coordinate along chromosome.
###     2) allelic_ratio_1 : allelic ratio of coordinate in project1 (child).
###     3) allelic_ratio_2 : allelic ratio of coordinate in project2 (parent).
###
### Uses genome definition files to only output data lines for chromosomes of interest.
### inputfile1 : 'SNP_CNV_v1.txt' file for experimental strain.
### inputFile2 : 'SNP_CNV_v1.txt' file for parent strain.
###

import string, sys, time

genome       = sys.argv[ 1]
genomeUser   = sys.argv[ 2]
project1     = sys.argv[ 3]
project1User = sys.argv[ 4]
project2     = sys.argv[ 5]
project2User = sys.argv[ 6]
main_dir     = sys.argv[ 7]

logName     = main_dir+"users/"+project1User+"/projects/"+project1+"/process_log.txt"
inputFile1  = main_dir+"users/"+project1User+"/projects/"+project1+"/SNP_CNV_v1.txt"
inputFile2  = main_dir+"users/"+project2User+"/projects/"+project2+"/SNP_CNV_v1.txt"
with open(logName, "a") as myfile:
	myfile.write("|\t    Comparing allelic ratios from '"+project1+"' and '"+project2+"',\n")

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("*====================================================================*\n")
	myfile.write("| Log of 'py/dataset_process_for_allelic_ratio_analysis.ddRADseq.py' |\n")
	myfile.write("*--------------------------------------------------------------------*\n")


#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+genomeUser+"/genomes/"+genome+"/"

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
	myfile.write("|\tProcessing parent and child  'putative_SNPs_v4.txt' file -> allelic fractions.\n")

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
	myfile.write("|\tGathering ratio data for each 'putative_SNPs_v4.txt' location.\n")

# Open dataset files.
data1 = open(inputFile1,"r")
data2 = open(inputFile2,"r")

with open(logName, "a") as myfile:
	myfile.write("|\t1")

print '### Data lines for each locus found in either the project or parent "SNP_CNV_v1.txt" files.'
print '###    chromosome_name'
print '###    bp_coordinate'
print '###    project1 allelic ratio: NA if not found, major bp_count divided by total count otherwise.'
print '###    project2 allelic ratio: NA if not found, major bp_count divided by total count otherwise.'
with open(logName, "a") as myfile:
	myfile.write("|\t2")

line1 = data1.readline()
line2 = data2.readline()
error_endOfFile = False
while (error_endOfFile == False):
	P1_chrID,P1_chrName,P1_position,P1_countA,P1_countT,P1_countG,P1_countC = processLine_CNV_SNP(line1)
	P2_chrID,P2_chrName,P2_position,P2_countA,P2_countT,P2_countG,P2_countC = processLine_CNV_SNP(line2)
	## Ensure that we're at the same chromosome in each file.
	#  Files are sorted by chromosome, with increasing chromosome IDs.
	#  If file load line hits the end of the file, then break out of the loops, as we're
	#      done with file processing.
	while (P1_chrID <> P2_chrID):
		while (P1_chrID < P2_chrID):
			# load a line from inputFile1.
			line1 = data1.readline()
			if not line1: # EOF 1
				error_endOfFile = True
				break
			P1_chrID,P1_chrName,P1_position,P1_countA,P1_countT,P1_countG,P1_countC = processLine_CNV_SNP(line1)
		while (P2_chrID < P1_chrID):
			# load a line from inputFile2.
			line2 = data2.readline()
			if not line2: # EOF 2
				loop_eror = True
				break
			P2_chrID,P2_chrName,P2_position,P2_countA,P2_countT,P2_countG,P2_countC = processLine_CNV_SNP(line2)
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
	while (P1_position <> P2_position):
		while (P1_position < P2_position):
			# load a line from inputFile1.
			line1 = data1.readline()
			if not line1: # endOfFile 1
				error_endOfFile = True
				break;
			P1_chrID,P1_chrName,P1_position,P1_countA,P1_countT,P1_countG,P1_countC = processLine_CNV_SNP(line1)
			if (P1_chrID <> P2_chrID): # endOfChr 1
				error_endOfChr = True
				break;
		while (P2_position < P1_position):
			# load a line from inputFile2.
			line2 = data2.readline()
			if not line2: # endOfFile 2
				loop_eror = True
				break
			P2_chrID,P2_chrName,P2_position,P2_countA,P2_countT,P2_countG,P2_countC = processLine_CNV_SNP(line2)
			if (P1_chrID <> P2_chrID): # endOfChr 2
				error_endOfChr = True
				break
		if error_endOfFile:
			break
		if error_endOfChr:
			break
	# End of while loop that registers file pointer to coordinate along chromosome.
	if error_endOfFile:
		break

	# At this time, line and line2 shuold point to the same chromosome and coordinate, or the
	# file processing loop will have exited because the locus was not found in one or the other
	# input data file.
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

	# Output new pileup line
	#     0) chr_name        : Chromosome name string.
	#     1) bp_coordinate   : bp coordinate along chromosome.
	#     2) allelic_ratio_1 : allelic ratio of coordinate in project1 (child).
	#     3) allelic_ratio_2 : allelic ratio of coordinate in project2 (parent).
	if P1_chrID!=0:
		print P1_chrName + '\t' + str(P1_position) + '\t' + str(P1_allelicRatio) + '\t' + str(P2_allelicRatio)
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
	myfile.write("*--------------------------------------------------------------------*\n")
	myfile.write("| End of 'py/dataset_process_for_allelic_ratio_analysis.ddRADseq.py' |\n")
	myfile.write("*====================================================================*\n")
