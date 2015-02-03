# Input arguments:
#	1) user    : Name of user account installing the genome file.
#	2) genome  : Name of the genome.
#	3) logfile : Path and file name of output log file.
#
# Output:


#------------------------------------------------------------------------------------------------------------
def find_nmer(seq):
	# Calculates position in nmer count vector and generates an error condition if the test sequence contains non-ATCG characters.
	steps = len(seq);
	nt    = [0]*steps;
	err   = 'false';
	seq   = seq.upper();
	for index in range(steps):
		if seq[index] == 'A':
			nt[index] = 0;
		elif seq[index] == 'T':
			nt[index] = 1;
		elif seq[index] == 'C':
			nt[index] = 2;
		elif seq[index] == 'G':
			nt[index] = 3;
		else:
			err = 'true';
	nmer = 1;
	for index in range(steps):
		nmer += nt[index]*4**index;
	return [nmer,err];
#------------------------------------------------------------------------------------------------------------
def rev_com(seq):
	# Generate the reverse complement sequence of an input DNA sequence string.
	seq         = seq.upper();
	rev_seq     = reversed(seq);
	rev_com_seq = rev_seq;
	for index in range(len(seq)):
		if   rev_com_seq[index] == 'A':
			rev_com_seq[index] = 'T';
		elif rev_com_seq[index] == 'T':
			rev_com_seq[index] = 'A';
		elif rev_com_seq[index] == 'C':
			rev_com_seq[index] = 'G';
		elif rev_com_seq[index] == 'G':
			rev_com_seq[index] = 'C';
		else:
			rev_com_seq[index] = 'n';
	return rev_com_seq;
#------------------------------------------------------------------------------------------------------------


import string, sys, re, time;
userName    = sys.argv[1];
genomeName  = sys.argv[2];
logName     = sys.argv[3];

# Initialize time counter and log file section.
t0 = time.clock();
with open(logName, "a") as myfile:
	myfile.write("\n\t\t-----------------------------------------------------------");
	myfile.write("\n\t\t Log of repetitiveness_1.py");
	myfile.write("\n\t\t-----------------------------------------------------------");


#============================================================================================================
# Determine name of reformatted reference FASTA file.
#------------------------------------------------------------------------------------------------------------

# Find name of genome FASTA file for species being examined.
#     Read in and parse : "Ymap_root/users/[userName]/genomes/[genomeName]/[genome]/reference.txt"
reference_file = workingDir + '/reference.txt';
refFile        = open(reference_file,'r');
refFASTA       = refFile.read().strip();
refFile.close();

# Generate the name of the reformatted FASTA file: "test.fasta" => "test.2.fasta"
refFASTA       = refFASTA.replace("fasta","2.fasta");
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tFASTA reformatted to single-line entries : " + refFASTA);


#============================================================================================================
# Determine which chromosomes from reference FASTA are being examined.
#------------------------------------------------------------------------------------------------------------

# Look up chromosome name strings for genome in use.
#     Read in and parse : "Ymap_root/users/[userName]/genomes/[genomeName]/figure_definitions.txt"
figureDefinition_file  = workingDir + '/figure_definitions.txt';
figureDefinitionFile   = open(figureDefinition_file,'r');
figureDefinitionData   = figureDefinitionFile.readlines();
#		Example lines in figureDefinition_file:
#			Chr  Use   Label   Name                         posX   posY   width   height
#			1    1     Chr1    Ca21chr1_C_albicans_SC5314   0.15   0.8    0.8     0.0625
#			2    1     Chr2    Ca21chr2_C_albicans_SC5314   0.15   0.7    *       0.0625
#			0    0     Mito    Ca19-mtDNA                   0.0    0.0    0.0     0.0

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tDetermining number of chromosomes of interest in genome.");

# Determine the number of chromosomes of interest in genome.
chrName_maxcount = 0;
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line));
	chr_num = line_parts[0];
	if chr_num.isdigit():
		chr_num    = int(float(line_parts[0]));
		chr_use    = int(float(line_parts[1]));
		chr_label  = line_parts[2];
		chr_name   = line_parts[3];
		if chr_num > 0:
			if chr_num > chrName_maxcount:
				chrName_maxcount = chr_num;
figureDefinitionFile.close();

# Pre-allocate chrName_array
chrName = [];
for x in range(0, chrName_maxcount):
	chrName.append([]);

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tGathering name strings for chromosomes.");

# Gather name strings for chromosomes, in order.
figureDefinitionFile  = open(figureDefinition_file,'r');
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line));
	chr_num = line_parts[0];
	if chr_num.isdigit():
		chr_num    = int(float(line_parts[0]));
		chr_use    = int(float(line_parts[1]));
		chr_label  = line_parts[2];
		chr_name   = line_parts[3];
		if chr_num != 0:
			chrName[int(float(chr_num))-1] = chr_name;
			with open(logName, "a") as myfile:
				myfile.write("\n\t\t\t\t\tChr" + str(chr_num) + " = " + chr_name);
figureDefinitionFile.close();

# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount;


#============================================================================================================
# Process used chromosomes into k-mer counts array.
#------------------------------------------------------------------------------------------------------------

# Initialize nmer_counts vector, zeros of length 4^nmer_length.
nmer_length = 10;
nmer_counts = [0]*(4**nmer_length);

# Open reformatted FASTA file.
FASTA_file = workingDir + refFASTA;
FASTA_data = open(FASTA_file,'r');

# Process reformatted FASTA file, entry by entry, to collect nmer counts.
while True:
	# FASTA entries are pairs of lines with the following structure.
	#    >Ca_a.chr1 (9638..10115) (478bp) [*]
	#    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
	line1 = FASTA_data.readline();
	line2 = FASTA_data.readline();
	if not line2:
		break  # EOF
	first_char = line1[:1];
	if first_char == ">":
		# First line is header to FASTA entry, so file contents are formatted properly.
		# Determine chromosome/contig name by isolating the first space-delimited string, then removing the header character ">".
		line_parts = string.split(string.strip(line1));
		chr_name   = line_parts[0];
		chr_name   = chr_name.replace(">","");

		# Run along chromosome, adding to k-mer count vector "nmer_counts".
		string_start = 0;
		string_end   = len(line2)-nmer_length;  # for a string of length 5, a nmer_length of 3 will have an end coordinate of 2 in [0, 1, 2].
		for index in range(string_end+1):
			test_string      = line2[index:(index+nmer_length)];
			results_1        = find_nmer(test_string);
			results_2        = find_nmer(rev_com(test_string));
			forward_nmer_num = results_1[0];
			forward_nmer_err = results_1[1];
			reverse_nmer_num = results_2[0];
			# If string is a valid DNA sequence:
			if forward_nmer_err == 'false':
				nmer_counts[forward_nmer_num] += 1;
				nmer_counts[reverse_nmer_num] += 1;


#============================================================================================================
# Process k-mer array into repetitiveness scores for each coordinate.
#------------------------------------------------------------------------------------------------------------
# chrName[] contains names of used chromosomes.
while True:
	# FASTA entries are pairs of lines with the following structure.
	#    >Ca_a.chr1 (9638..10115) (478bp) [*]
	#    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
	line1 = FASTA_data.readline();
	line2 = FASTA_data.readline();
	if not line2:
		break  # EOF
	first_char = line1[:1];
	if first_char == ">":
		# First line is header to FASTA entry, so file contents are formatted properly.
		# Determine chromosome/contig name by isolating the first space-delimited string, then removing the header character ">".
		line_parts = string.split(string.strip(line1));
		chr_name   = line_parts[0];
		chr_name   = chr_name.replace(">","");

		# Run along chromosome, adding to k-mer count vector "nmer_counts".
		string_start = 0;
		string_end   = len(line2)-nmer_length;  # for a string of length 5, a nmer_length of 3 will have an end coordinate of 2 in [0, 1, 2].
		for index in range(string_end+1):
			test_string      = line2[index:(index+nmer_length)];
			results_1        = find_nmer(test_string);
			results_2        = find_nmer(rev_com(test_string));
			forward_nmer_num = results_1[0];
			forward_nmer_err = results_1[1];
			reverse_nmer_num = results_2[0];
			# If string is a valid DNA sequence:
			if forward_nmer_err == 'false':
				nmer_counts[forward_nmer_num] += 1;
				nmer_counts[reverse_nmer_num] += 1;


#============================================================================================================
# 
#------------------------------------------------------------------------------------------------------------


#============================================================================================================
# 
#------------------------------------------------------------------------------------------------------------


#============================================================================================================
# Code section to output information about genome restriction fragments.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tOutputting repetitiveness data of standard-bin fragmented genome.")

for fragment in range(1,numFragments):
	# Output a line for each fragment.
	#     fragments[fragment-1] = [chr_num,bp_start,bp_end, data_count,data_max,length]
	#     0) chr_num
	#     1) bp_start
	#     2) bp_end
	#     3) repetitiveness_max
	#     4) repetitiveness_ave
	#     5) length = bp_end-bp_start+1

	chr_num         = fragments[fragment-1][0]
	bp_start        = fragments[fragment-1][1]
	bp_end          = fragments[fragment-1][2]

	repet_max       = fragments[fragment-1][4]
	repet_ave       = fragments[fragment-1][5]
	fragment_length = bp_end - bp_start + 1

	print str(chr_num) + '\t' + str(bp_start) + '\t' + str(bp_end) + '\t' + str(repet_max) + '\t' + str(repet_ave) + '\t' + str(fragment_length)

#------------------------------------------------------------------------------------------------------------
# End of code section to output information about fragments. 
#============================================================================================================
