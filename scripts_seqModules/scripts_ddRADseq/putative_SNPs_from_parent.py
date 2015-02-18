###
### Simplify parental putative_SNP list to contain only those loci with an allelic ratio on range [0.25 .. 0.75].
###
### Uses genome definition files to only output data lines for chromosomes of interest.
###


def process_ParentLine(entry_line):
	global chrNums
	global chrName
	global chrCount
	# Process 'putative_SNPs_v4.txt' file line.
	# example lines:
	#       chromosome                   coord   ref   A    T     G   C
	#       Ca21chr1_C_albicans_SC5314   13988   T     1    12    0   0
	#       Ca21chr1_C_albicans_SC5314   13993   A     12   0     0   1
	#       Ca21chr1_C_albicans_SC5314   14003   T     1    412   0   0
	#       Ca21chr1_C_albicans_SC5314   14004   T     1    413   0   0
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

import string, sys, time

genome            = sys.argv[ 1]
genomeUser        = sys.argv[ 2]
projectChild      = sys.argv[ 3]
projectChildUser  = sys.argv[ 4]
projectParent     = sys.argv[ 5]
projectParentUser = sys.argv[ 6]
main_dir          = sys.argv[ 7]

logName           = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/process_log.txt"
inputFile_P       = main_dir+"users/"+projectParentUser+"/projects/"+projectParent+"/putative_SNPs_v4.txt"
inputFile_C       = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/SNP_CNV_v1.txt"

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\t\t*===================================================*\n")
	myfile.write("\t\t| Log of 'scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_parent.py'          |\n")
	myfile.write("\t\t*---------------------------------------------------*\n")


#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+genomeUser+"/genomes/"+genome+"/"

#============================================================================================================
# Load FastaName from 'reference.txt' for genome in use.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t|\tIdentifying name of reference FASTA file.\n")
reference_file = genomeDirectory + '/reference.txt'
refFile        = open(reference_file,'r')
FastaName      = refFile.read().strip()
refFile.close()
FastaName      = FastaName.replace(".fasta", "")

#============================================================================================================
# Process 'preprocessed_SNPs.txt' file for projectParent to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t|\tProcessing parent 'putative_SNPs_v4' file -> het loci.\n")

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
	myfile.write("\t\t|\tDetermining number of chromosomes of interest in genome.\n")

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
	myfile.write("\t\t|\tGathering name strings for chromosomes.\n")

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
				myfile.write("\t\t|\t\t" + str(chr_num) + " : " + chr_name + " = " + chr_nameShort + "\n")
figureDefinitionFile.close()

# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount
with open(logName, "a") as myfile:
	myfile.write("\t\t|\t\tMax chr string : "+str(chrCount)+"\n")
#............................................................................................................

with open(logName, "a") as myfile:
	myfile.write("\t\t|\tOpen parent 'putative_SNPs_v4.txt' file.\n")

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
	myfile.write("\t\t|\tGathering read coverage data for each fragment.\n")

# Open dataset 'putative_CNVs_v1.txt' file.
data_P = open(inputFile_P,"r")

print '### Data lines for each het locus in parent : [chromosome_name, bp_coordinate, countA, countT, countG, countC]'

# Process 'SNP_CNV_v1.txt' file for both parents, line by line... while checking for missing data.
line_P = data_P.readline()
error_endOfFile = False
while (error_endOfFile == False):
	P_chrID,P_chrName,P_position,P_countA,P_countT,P_countG,P_countC = process_ParentLine(line_P)
	P_list = [int(float(P_countA)), int(float(P_countT)), int(float(P_countG)), int(float(P_countC))]
	if (sum(P_list) == 0):
		P_allelicRatio = 0
	else:
		P_allelicRatio = max(P_list)/float(sum(P_list))
	if ((P_allelicRatio < 0.75) and (P_allelicRatio > 0.25) and (sum(P_list) > 20)):
		print P_chrName+"\t"+str(P_position)+"\t"+str(P_countA)+"\t"+str(P_countT)+"\t"+str(P_countG)+"\t"+str(P_countC)
	error_endOfFile = False
	line_P = data_P.readline()
	if not line_P: # EOF 1
		error_endOfFile = True
		break
data_P.close()

#------------------------------------------------------------------------------------------------------------
# End of main code block.
#============================================================================================================

print '### End of preprocessed parental SNP, child SNP data.'

with open(logName, "a") as myfile:
	myfile.write("\t\t|\tTime to process = " + str(time.clock()-t0) + "\n")
	myfile.write("\t\t*---------------------------------------------------*\n")
	myfile.write("\t\t| End of 'py/putative_SNPs_from_parent.py'          |\n")
	myfile.write("\t\t*===================================================*\n")
