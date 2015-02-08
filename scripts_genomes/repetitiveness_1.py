# Input arguments:
#	1) user     : Name of user account installing the genome file.
#	2) genome   : Name of the genome.
#	3) main_dir : Root directory.
#	4) logfile  : Path and file name of output log file.
#
import string, sys, re, time;
userName    = sys.argv[1];
genomeName  = sys.argv[2];
main_dir    = sys.argv[3];
logName     = sys.argv[4];

nmer_length = 10;

with open(logName, "a") as myfile:
	myfile.write("\t\t\t*================================================================*\n")
	myfile.write("\t\t\t| Log of 'repetitiveness_1.py'                                   |\n")
	myfile.write("\t\t\t*----------------------------------------------------------------*\n")


#------------------------------------------------------------------------------------------------------------
# Calculates position in nmer count vector and generates an error condition if the test sequence contains non-ATCG characters.
def find_nmer(seq):
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
	nmer = 0;  # vectors in python are 0-rooted.
	for index in range(steps):
		nmer += nt[index]*4**index;
	return [nmer,err];

#------------------------------------------------------------------------------------------------------------
# Recursively reverse the input string.
def reverse(text):
	if len(text) <= 1:
		return text;
	return reverse(text[1:]) + text[0];

#------------------------------------------------------------------------------------------------------------
# Generate the reverse complement sequence of an input DNA sequence string.
def rev_com(seq):
	seq         = seq.upper();
	rev_seq     = reverse(seq);
	rev_com_seq = '';
	for index in range(len(seq)):
		if   rev_seq[index] == 'A':
			rev_com_seq += 'T';
		elif rev_seq[index] == 'T':
			rev_com_seq += 'A';
		elif rev_seq[index] == 'C':
			rev_com_seq += 'G';
		elif rev_seq[index] == 'G':
			rev_com_seq += 'C';
		else:
			rev_com_seq += 'n';
	return rev_com_seq;

#------------------------------------------------------------------------------------------------------------


# Initialize time counter and log file section.
t0 = time.clock();


#============================================================================================================
# Determine name of reformatted reference FASTA file.
#------------------------------------------------------------------------------------------------------------
# Find name of genome FASTA file for species being examined.
#     Read in and parse : "Ymap_root/users/[userName]/genomes/[genomeName]/[genome]/reference.txt"
workingDir     = main_dir + 'users/' + userName + '/genomes/' + genomeName + '/';
reference_file = workingDir + 'reference.txt';
refFile        = open(reference_file,'r');
refFASTA       = refFile.read().strip();
refFile.close();

# Generate the name of the reformatted FASTA file: "test.fasta" => "test.2.fasta"
refFASTA       = refFASTA.replace(".fasta",".2.fasta");
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
#	Example lines in figureDefinition_file:
#		Chr  Use   Label   Name                         posX   posY   width   height
#		1    1     Chr1    Ca21chr1_C_albicans_SC5314   0.15   0.8    0.8     0.0625
#		2    1     Chr2    Ca21chr2_C_albicans_SC5314   0.15   0.7    *       0.0625
#		0    0     Mito    Ca19-mtDNA                   0.0    0.0    0.0     0.0

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
nmer_counts = [0]*(4**nmer_length);

# Open reformatted FASTA file.
FASTA_file = workingDir + refFASTA;
FASTA_data = open(FASTA_file,'r');

# Reset timer.
t0a = time.clock();
t1  = t0a;
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t#### Start tallying nmers.\n");

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
	old_chr_name = '';
	if first_char == ">":
		# First line is header to FASTA entry, so file contents are formatted properly.
		# Determine chromosome/contig name by isolating the first space-delimited string, then removing the header character ">".
		line_parts = string.split(string.strip(line1));
		chr_name   = line_parts[0];
		chr_name   = chr_name.replace(">","");

		# If the current chromosome is different than the last one, output log entry.
		if (chr_name <> old_chr_name):
			with open(logName, "a") as myfile:
				myfile.write("\t\t\t####\t" + str(time.clock() - t1) + " seconds for this chr.\n");
				t1 = time.clock();
				if chr_name in chrName:
					myfile.write("\t\t\t#### Tallying nmers of: " + chr_name + "\n");
				else:
					myfile.write("\t\t\t#### Skipping '" + chr_name + "' because it was deselected.\n");

		# If the current chromosome is one of those in use...
		if chr_name in chrName:
			# Run along chromosome coordinates, from start(0) to end (len(line2)-nmer_length+1), such that a nmer can be examined starting at each coordinate.
			# For each coordinate, add to k-mer count vector "nmer_counts".
			count = 0;
			for index in range(0, len(line2)-nmer_length+1):
				count += 1;
				# current nmer string.
				test_string      = line2[index:(index+nmer_length)];
				# determine nmer_counts vector position for the test_string and its reverse complement.
				results_1        = find_nmer(test_string);
				results_2        = find_nmer(rev_com(test_string));
				forward_nmer_num = results_1[0];
				forward_nmer_err = results_1[1];
				reverse_nmer_num = results_2[0];
				if forward_nmer_err == 'false':
					# If test_string is a valid DNA sequence, increment the nmer_counts vector position for the string and its reverse complement.
					nmer_counts[forward_nmer_num] += 1;
					nmer_counts[reverse_nmer_num] += 1;
		old_chr_name = chr_name;
FASTA_data.close();

# Reset timer.
with open(logName, "a") as myfile:
	if old_chr_name in chrName:
		myfile.write("\t\t\t####\t" + str(time.clock() - t1) + " seconds for this chr.\n");
	myfile.write("\t\t\t####\n");
	myfile.write("\t\t\t#### " + str(time.clock() - t0a) + " seconds to talley incidence of nmers across genome.\n");
	myfile.write("\t\t\t####\n");
t0b = time.clock();
t1  = t0b;
with open(logName, "a") as myfile:
	myfile.write("\t\t\t#### Start generating repetitiveness score across genome.\n");


#============================================================================================================
# Output header lines for repetititveness file output.
#------------------------------------------------------------------------------------------------------------
print '## Repetitiveness score per bp location.';
print '##';
print '## columns = [chrName, bpCoordinate, repetitivenessScore]';


#============================================================================================================
# Process k-mer array into repetitiveness scores for each coordinate, then generate output lines.
#------------------------------------------------------------------------------------------------------------
FASTA_data = open(FASTA_file,'r');
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
	old_chr_name = '';
	if first_char == ">":
		# First line is header to FASTA entry, so file contents are formatted properly.
		# Determine chromosome/contig name by isolating the first space-delimited string, then removing the header character ">".
		line_parts = string.split(string.strip(line1));
		chr_name   = line_parts[0];
		chr_name   = chr_name.replace(">","");

		# If the current chromosome is different than the last one, output log entry.
		if (chr_name <> old_chr_name):
			with open(logName, "a") as myfile:
				myfile.write("\t\t\t####\t" + str(time.clock() - t1) + " seconds for this chr.\n");
			t1 = time.clock();
			with open(logName, "a") as myfile:
				myfile.write("\t\t\t#### Outputing repetitiveness scores for: " + chr_name + "\n");

		# If the current chromosome is one of those in use...
		if chr_name in chrName:
			# Run along chromosome coordinates, from start (0) to end (len(line2)-1), such that all nmers overlapping the coordinate can be examined.
			for index in range(0, len(line2)-1):
				score_sum = 0;
				for index_offset in range(-nmer_length+1,1):
					# current nmer string.
					if (index+index_offset >= 0) and (index+index_offset+nmer_length-1 <= len(line2)-1):
						test_string      = line2[(index+index_offset):(index+index_offset+nmer_length)];
						# determine nmer_counts vector position for the test_string. The reverse-complement was dealt with in the tallying step previous.
						results_1        = find_nmer(test_string);
						forward_nmer_num = results_1[0];
						forward_nmer_err = results_1[1];
						if forward_nmer_err == 'false':
							# If test_string is a valid DNA sequence, add to score_sum.
							score_sum += nmer_counts[forward_nmer_num];
				# Output repetitiveness score line to file: [chrName, bpCoordinate, repetitivenessScore]
				print chr_name + '\t' + str(index+1) + '\t' + str(score_sum);
		old_chr_name = chr_name;


#============================================================================================================
# Conclude log outputs.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t\t####\t" + str(time.clock() - t1) + " seconds for this chr.\n");
	myfile.write("\t\t\t####\n");
	myfile.write("\t\t\t#### Total time for computation of genome repetitiveness = " + str(time.clock() - t0) + "\n");
	myfile.write("\t\t\t####\n");
