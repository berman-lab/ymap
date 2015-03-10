### 
### Simplify child putative_SNP list to contain only those loci found in the hapmap.
###
### Uses genome definition files to only output data lines for chromosomes of interest.
###


def process_ChildLine(entry_line):
	global chrNums;
	global chrName;
	global chrCount;
	# Process 'SNP_CNV_v1.txt' file line.
	# example lines:
	#       chromosome               coord   total   ref   A    T    G    C
	#       ChrA_C_glabrata_CBS138   45      36      C     0    0    0    36
	#       ChrA_C_glabrata_CBS138   46      37      T     0    37   0    0
	#       ChrA_C_glabrata_CBS138   47      38      A     38   0    0    0
	#       ChrA_C_glabrata_CBS138   48      39      A     39   0    0    0
	child_line = string.strip(entry_line);
	child_line = child_line.split('\t');
	C_chr_name = child_line[0];   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	C_position = child_line[1];   # chr position of bp.      : 2286371
	C_countTot = child_line[2];   # total count at bp.       : 101
	C_refBase  = child_line[3];   # reference base at bp.    : T
	C_countA   = child_line[4];   # count of A.              : 100
	C_countT   = child_line[5];   # count of T.              : 0
	C_countG   = child_line[6];   # count of G.              : 0
	C_countC   = child_line[7];   # count of C.              : 1
	# Determine chrID associated with chromosome name.
	C_chr = 0;
	for x in range(0,chrCount):
		if (chrNums[x] != 0):
			if chrName[x] == C_chr_name:
				C_chr = x+1;
	C_chrName = chrName[C_chr-1];
	return C_chr,C_chrName,C_position,C_countA,C_countT,C_countG,C_countC;

def process_HapmapLine(entry_line):
	global chrNums;
	global chrName;
	global chrCount;
	# Process 'SNPdata_parent.txt' file line.
	# example lines:
	#       chromosome                   coord   HomA   HomB   Status1   (Status2 ...)
	#       Ca21chr1_C_albicans_SC5314   812     C      T      0         (1       ...)
	#       Ca21chr1_C_albicans_SC5314   816     T      C      0         (1       ...)
	#       Ca21chr1_C_albicans_SC5314   879     G      A      0         (0       ...)
	hapmap_line   = string.strip(entry_line);
	hapmap_line   = hapmap_line.split('\t');
	H_chr_name    = hapmap_line[0];   # chromosome   : Ca21chrR_C_albicans_SC5314
	H_position    = hapmap_line[1];   # coordinate   : 2286371
	H_status_list = [];
	for entry in range(4,(len(hapmap_line))):
		H_status_list.append(int(hapmap_line[entry]));
		# entry values
		#    0  : correct phase.
		#    1  : incorrect phase.
		#    10 : no phase information, due to no data at coodinate in child dataset.
		#    11 : no phase information, due to coordinate not matching a LOH fragment in child dataset.
		#    12 : no phase information, due to surprise allele in child dataset.

	## Determine consensus hapmap entry.
	if (H_status_list.count(1) == 0) and (H_status_list.count(0) == 0):
		# No useful hapmap entries for this locus...
		if (H_status_list.count(11) > H_status_list.count(10)) and (H_status_list.count(11) > H_status_list.count(12)):
			# ...because the coordinate was not in defined LOH regions during setup of hapmap.
			H_status = 11;
		elif (H_status_list.count(10) > H_status_list.count(11)) and (H_status_list.count(10) > H_status_list.count(12)):
			# ...because the coordinate had no data during setup of hapmap.
			H_status = 10;
		elif (H_status_list.count(12) > H_status_list.count(10)) and (H_status_list.count(12) > H_status_list.count(11)):
			# ...because the coordinate had a surprise allele during setup of hapmap.
			H_status = 12;
		else:
			# ...because some undetermined combination of errors at this coordinate.
			H_status = 13;
	elif (H_status_list.count(1) > H_status_list.count(0)):
		# More hapmap entries indicate phasing 1 for this locus.
		H_status = 1;
	elif (H_status_list.count(1) == H_status_list.count(0)):
		# Useful hapmap entries are ambiguous, so random phasing is chosen.
		H_status = random.randrange(0,2);
	else:
		# More hapmap entries indicate phasing 0 for this locus.
		H_status = 0;
	# Determine chrID associated with chromosome name.
	H_chr = 0;
	for x in range(0,chrCount):
		if (chrNums[x] != 0):
			if chrName[x] == H_chr_name:
				H_chr = x+1;
	H_chrName = chrName[H_chr-1];
	return H_chr,H_chrName,H_position,H_status;

import string, sys, time, random;
random.seed();

genome             = sys.argv[1];
genomeUser         = sys.argv[2];
projectChild       = sys.argv[3];
projectChildUser   = sys.argv[4];
hapmap             = sys.argv[5];
HapmapUser         = sys.argv[6];
main_dir           = sys.argv[7];

logName            = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/process_log.txt";
inputFile_H        = main_dir+"users/"+HapmapUser+"/hapmaps/"+hapmap+"/SNPdata_parent.txt";
inputFile_C        = main_dir+"users/"+projectChildUser+"/projects/"+projectChild+"/SNP_CNV_v1.txt";

t0 = time.clock();

with open(logName, "a") as myfile:
	myfile.write("\t\t*====================================================================================*\n");
	myfile.write("\t\t| Log of 'scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_hapmap_in_child.py' |\n");
	myfile.write("\t\t*------------------------------------------------------------------------------------*\n");


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
# Process 'preprocessed_SNPs.txt' file for hapmap to determine initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t|\tProcessing hapmap file -> het loci.\n")

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


# Process hapmap file, as well as "SNP_CNV_v1.txt" for the data from the child.
with open(logName, "a") as myfile:
        myfile.write("\t\t|\tLoading SNP coordinates from hapmap.\n");
print '### Data lines for each locus in hapmap : [chromosome_name, bp_coordinate, countA, countT, countG, countC]';
data_H               = open(inputFile_H,"r");
old_H_chrID          = 0;
hapmap_loci          = [];
hapmap_loci_unphased = [];
for line_H in data_H:
	if (len(line_H) > 0):
		if (line_H[0] != "#"):
			H_chrID,H_chrName,H_position,H_status = process_HapmapLine(line_H);
			if H_chrID != old_H_chrID:
				with open(logName, "a") as myfile:
					myfile.write("\t\t|\t\tchr = "+str(H_chrName)+"\n");
			if H_chrID > 0:
				if H_status in [0, 1]:
					hapmap_loci.append([H_chrName,H_position]);
				elif H_status in [10,11,12,13]:
					hapmap_loci_unphased.append([H_chrName,H_position]);
			old_H_chrID = H_chrID;
data_H.close();

# Process child dataset for matches to hapmap.
with open(logName, "a") as myfile:
	myfile.write("\t\t|\tScreening through child dataset for coordinates matching hapmap loci.");
data_C                    = open(inputFile_C,"r");
old_C_chrID               = 0;
child_SNPs                = [];
child_SNPs_small          = [];
child_SNPs_unphased       = [];
child_SNPs_unphased_small = [];
counter                   = 0;
for line_C in data_C:
	if (len(line_C) > 0):
		if (line_C[0] != "#"):
			C_chrID,C_chrName,C_position,C_countA,C_countT,C_countG,C_countC = process_ChildLine(line_C)
			if C_chrID != old_C_chrID:
				with open(logName, "a") as myfile:
					myfile.write("\n\t\t|\t\tchr = "+str(C_chrName)+"\n");
				counter = 0;
			if (C_chrID > 0) and (int(C_countA)+int(C_countT)+int(C_countG)+int(C_countC) >= 20):   # chromosome is identified and in use; read depth >= 20.
				testval = 0;
				if [C_chrName,C_position] in hapmap_loci:
					child_SNPs.append([C_chrName,C_position,C_countA,C_countT,C_countG,C_countC]);
					child_SNPs_small.append([C_chrName,C_position]);
					testval = 1;
				if [C_chrName,C_position] in hapmap_loci_unphased:
					child_SNPs_unphased.append([C_chrName,C_position,C_countA,C_countT,C_countG,C_countC]);
					child_SNPs_unphased_small.append([C_chrName,C_position]);
					testval = 1;
				if (testval == 1):
					if counter == 0:
						with open(logName, "a") as myfile:
							myfile.write("\t\t|\t\t");
					if counter%10 == 0:
						with open(logName, "a") as myfile:
							myfile.write(".");
					if counter == 800:
						with open(logName, "a") as myfile:
							myfile.write("\n\t\t|\t\t");
						counter = 0;
					counter += 1;
			old_C_chrID = C_chrID;
data_C.close();

# Output child lines from hapmap positions.
with open(logName, "a") as myfile:
	myfile.write("\n\t\t|\tOutputting lines from child dataset that match coordinates of hapmap loci.\n");
for SNP in hapmap_loci:
	if SNP in child_SNPs_small:
		SNP_data = child_SNPs[child_SNPs_small.index(SNP)];
		print SNP_data[0]+"\t"+SNP_data[1]+"\t"+SNP_data[2]+"\t"+SNP_data[3]+"\t"+SNP_data[4]+"\t"+SNP_data[5]+"\t+";
	if SNP in child_SNPs_unphased_small:
		SNP_data = child_SNPs_unphased[child_SNPs_unphased_small.index(SNP)];
		print SNP_data[0]+"\t"+SNP_data[1]+"\t"+SNP_data[2]+"\t"+SNP_data[3]+"\t"+SNP_data[4]+"\t"+SNP_data[5]+"\t-";

#------------------------------------------------------------------------------------------------------------
# End of main code block.
#============================================================================================================

print '### End of preprocessed hapmap loci vs. child SNP data.'

with open(logName, "a") as myfile:
	myfile.write("\t\t|\tTime to process = " + str(time.clock()-t0) + "\n")
	myfile.write("\t\t*------------------------------------------------------------------------------------*\n")
	myfile.write("\t\t| End of 'scripts_seqModules/scripts_ddRADseq/putative_SNPs_from_hapmap_in_child.py' |\n")
	myfile.write("\t\t*====================================================================================*\n")
