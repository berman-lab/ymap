### 
### Simplify child putative_SNP list to contain only those loci with an allelic ratio on range [0.25 .. 0.75] in the parent dataset.
###
### Uses genome definition files to only output data lines for chromosomes of interest.
###


def process_ChildLine(entry_line):
	global chrNums
	global chrName
	global chrCount
	# Process 'SNP_CNV_v1.txt' file line.
	# example lines:
	#       chromosome               coord   total   ref   A    T    G    C
	#       ChrA_C_glabrata_CBS138   45      36      C     0    0    0    36
	#       ChrA_C_glabrata_CBS138   46      37      T     0    37   0    0
	#       ChrA_C_glabrata_CBS138   47      38      A     38   0    0    0
	#       ChrA_C_glabrata_CBS138   48      39      A     39   0    0    0
	child_line = string.strip(entry_line)
	child_line = child_line.split('\t')
	C_chr_name = child_line[0]   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	C_position = child_line[1]   # chr position of bp.      : 2286371
	C_countTot = child_line[2]   # total count at bp.       : 101
	C_refBase  = child_line[3]   # reference base at bp.    : T
	C_countA   = child_line[4]   # count of A.              : 100
	C_countT   = child_line[5]   # count of T.              : 0
	C_countG   = child_line[6]   # count of G.              : 0
	C_countC   = child_line[7]   # count of C.              : 1
	# Determine chrID associated with chromosome name.
	C_chr = 0
	for x in range(0,chrCount):
		if (chrNums[x] != 0):
			if chrName[x] == C_chr_name:
				C_chr = x+1
	C_chrName = chrName[C_chr-1]
	return C_chr,C_chrName,C_position,C_countA,C_countT,C_countG,C_countC

def process_trimmedParentLine(entry_line):
	global chrNums
	global chrName
	global chrCount
	# Process 'trimmed_SNPs_v4.parent.txt' file line.
	# example lines:
	#       chromosome                   coord   A    T     G   C
	#       Ca21chr1_C_albicans_SC5314   13988   1    12    0   0
	#       Ca21chr1_C_albicans_SC5314   13993   12   0     0   1
	#       Ca21chr1_C_albicans_SC5314   14003   1    412   0   0
	#       Ca21chr1_C_albicans_SC5314   14004   1    413   0   0
	parent_line = string.strip(entry_line)
	parent_line = parent_line.split('\t')
	P_chr_name  = parent_line[0]   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	P_position  = parent_line[1]   # chr position of bp.      : 2286371
	# Determine chrID associated with chromosome name.
	P_chr = 0
	for x in range(0,chrCount):
		if (chrNums[x] != 0):
			if chrName[x] == P_chr_name:
				P_chr = x+1
	P_chrName = chrName[P_chr-1]
	return P_chr,P_chrName,P_position

import string, sys, time

genome             = sys.argv[1];
genomeUser         = sys.argv[2];
projectChild       = sys.argv[3];
projectChildUser   = sys.argv[4];
main_dir           = sys.argv[5];

logName            = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/process_log.txt";
inputFile_trimmedP = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/trimmed_SNPs_v4.parent.txt";
inputFile_C        = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/SNP_CNV_v1.txt";

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\t\t*======================================================================================*\n")
	myfile.write("\t\t| Log of 'scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_parent_in_child.2.py' |\n")
	myfile.write("\t\t*--------------------------------------------------------------------------------------*\n")


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
# Process 'preprocessed_SNPs.ddRADseq.txt' file for project to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t|\tProcessing 'putative_SNPs_v4' file -> het loci.\n")

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
	myfile.write("\t\t|\tOpen 'trimmed_SNPs_v4.parent.txt' file.\n")

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
	myfile.write("\t\t|\tCollecting data in child from putatitve SNP loci in parent.\n")

# Open dataset 'trimmed_SNPs_v4.parent.txt' file.
data_P = open(inputFile_trimmedP,"r")

print '### Data lines for each het locus in parent : [chromosome_name, bp_coordinate, countA, countT, countG, countC]'

# Process "trimmed_SNPs_v4.parent.txt" file containing SNP position data from the parent, as well as "SNP_CNV_v1.txt" for the data from the child.
line_P = data_P.readline();
error_endOfFile = False;
old_P_chrID = 0;
while (error_endOfFile == False):
	P_chrID,P_chrName,P_position = process_trimmedParentLine(line_P);
	if P_chrID != old_P_chrID:
		with open(logName, "a") as myfile:
			myfile.write("\t\t|\t\tchr = "+str(P_chrName)+"\n");

	data_C = open(inputFile_C,"r")
	line_C = data_C.readline()
	error_endOfFile = False
	while (error_endOfFile == False):
		C_chrID,C_chrName,C_position,C_countA,C_countT,C_countG,C_countC = process_ChildLine(line_C)
		if (C_chrID == P_chrID) and (C_position == P_position):
			print C_chrName+"\t"+str(C_position)+"\t"+str(C_countA)+"\t"+str(C_countT)+"\t"+str(C_countG)+"\t"+str(C_countC)
			data_C.close()
			break;
		line_C = data_C.readline()
		if not line_C: # EOF 2
			error_endOfFile = True
			data_C.close()
			break;
	if (error_endOfFile == True):
		# data line corresponding to parent SNP was not found in child "SNP_CNV_v1.txt" file.
		print P_chrName+"\t"+str(P_position)+"\t0\t0\t0\t0"
	data_C.close()

	error_endOfFile = False
	line_P = data_P.readline()
	if not line_P: # EOF 1
		error_endOfFile = True
		break
	old_P_chrID = P_chrID;
data_P.close()

#------------------------------------------------------------------------------------------------------------
# End of main code block.
#============================================================================================================

print '### End of preprocessed parental SNP, child SNP data.'

with open(logName, "a") as myfile:
	myfile.write("\t\t|\tTime to process = " + str(time.clock()-t0) + "\n")
	myfile.write("\t\t*--------------------------------------------------------------------------------------*\n")
	myfile.write("\t\t| End of 'scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_parent_in_child.2.py' |\n")
	myfile.write("\t\t*======================================================================================*\n")
