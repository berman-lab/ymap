# Input arguments:
#	1) user     : Name of user account installing the genome file.
#	2) genome   : Name of installed genome.
#	3) main_dir : Root directory.
#	4) logfile  : Path and file name of log file.
#
# Processing:
#	1) Load genome FASTA.
#	2) Output standard-bin fragmented genome FASTA.
#
import string, sys, re, time;
userName    = sys.argv[1];
genomeName  = sys.argv[2];
main_dir    = sys.argv[3];
logName     = sys.argv[4];


# Initialize time counter and log file section.
t0 = time.clock();
with open(logName, "a") as myfile:
	myfile.write("\n\t\t-----------------------------------------------------------");
	myfile.write("\n\t\t Log of genome_process_for_standard_bins_1.py");
	myfile.write("\n\t\t-----------------------------------------------------------");


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
#       Example lines in figureDefinition_file:
#               Chr  Use   Label   Name                         posX   posY   width   height
#               1    1     Chr1    Ca21chr1_C_albicans_SC5314   0.15   0.8    0.8     0.0625
#               2    1     Chr2    Ca21chr2_C_albicans_SC5314   0.15   0.7    *       0.0625
#               0    0     Mito    Ca19-mtDNA                   0.0    0.0    0.0     0.0

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
# Determine chromosome sizes and standard bin size.
#------------------------------------------------------------------------------------------------------------
# standard bin size is round(max(chr_lengths)/700).
#------------------------------------------------------------------------------------------------------------
# Open reformatted FASTA file.
FASTA_file = workingDir + refFASTA;
FASTA_data = open(FASTA_file,'r');

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t#### Determine chromosome lengths.");

# Process reformatted FASTA file, entry by entry, to determine maximum chromosome length.
chr_lenths = [];
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
		# If the current chromosome is one of those in use...
		if chr_name in chrName:
			chr_lenths.append(len(line2));
FASTA_data.close();
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t####     max_chr_length         = " + str(max(chr_lenths)) + " bp");

# Calculate standard bin size : round(max(chr_lenths)/700);
bases_per_bin = int(round(max(chr_lenths)/700));
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t####     bases_per_standard_bin = " + str(bases_per_bin) + " bp");

# Reset timer.
t0a = time.clock();
t1  = t0a;
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t#### Start genome fragmentation.\n");


#============================================================================================================
# Digest chromosome sequences and output new FASTA file entries.
#------------------------------------------------------------------------------------------------------------
# Open reformatted FASTA file.
FASTA_data = open(FASTA_file,'r');
# Process reformatted FASTA file, entry by entry, to digest into fragments.
chr_num = 0;
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
		# Increment chromosome number counter.
		chr_num += 1;

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
					myfile.write("\t\t\t#### Fragmenting: " + chr_name + "\n");
				else:
					myfile.write("\t\t\t#### Not fragmenting: " + chr_name + "\n");

		# If the current chromosome is one of those in use...
		if chr_name in chrName:
			fragment_start = 1;
			fragment_end   = bases_per_bin;
			while fragment_start < len(line2):
				fragment_string = line2[(fragment_start-1):(min(fragment_end,len(line2))-1)];
				fragment_length = min(fragment_end,len(line2)) - fragment_start + 1;
				print ">" + genomeName + ".chr" + str(chr_num) + " (" + str(fragment_start) + ".." + str(min(fragment_end,len(line2))) + ") (" + str(fragment_length) + "bp) [*]";
				print fragment_string;
				fragment_start += bases_per_bin;
				fragment_end   += bases_per_bin;
FASTA_data.close();
