# Input arguments: (Those with '[*]' at end are used here.)
# Process genome repetitiveness file by restriction fragment definitions.
#	1) test.repetitiveness.txt		(Ymap_root/users/[user]/genomes/[genome]/)
#	2) test.MfeI_MboI.fasta			(Ymap_root/users/[user]/genomes/[genome]/)
# Output fragment repetitiveness data.
#	1) test.repetitiveness.MfeI_MboI.txt	(Ymap_root/users/[user]/genomes/[genome]/)
#



import string, sys, re, time;
userName    = sys.argv[1];
genomeName  = sys.argv[2];
main_dir    = sys.argv[3];
logName     = sys.argv[4];

t0 = time.clock();
with open(logName, "a") as myfile:
	myfile.write("\t\t\t*================================================================*\n");
	myfile.write("\t\t\t| Log of 'genome_process_for_RNAseq.repetitiveness_2.py'         |\n");
	myfile.write("\t\t\t*----------------------------------------------------------------*\n");


#============================================================================================================
# Process restriction-digested genome file.
#------------------------------------------------------------------------------------------------------------
# Example FASTQ header line.
#     >Ca_a.chr1 (9638..10115) (478bp) [*]
# Lines ending in '[*]' are usable.   Other lines aren't usable.

# Find name of genome FASTA file for species being examined.
#     Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/reference.txt"
workingDir     = main_dir + 'users/' + userName + '/genomes/' + genomeName + '/';
reference_file = workingDir + '/reference.txt';
refFile        = open(reference_file,'r');
refFASTA       = refFile.read().strip();
refFile.close();

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tReference FASTA : " + refFASTA);

# Open RNAseq-digested genome FASTQ file.
expressionGenome    = refFASTA.replace(".fasta",".expression.fasta");
datafile_name       = refFASTA.replace(".fasta",".repetitiveness.txt");
datafile            = workingDir + datafile_name;
RNAseq_FASTA_file   = workingDir + expressionGenome;
RNAseq_FASTA_data   = open(RNAseq_FASTA_file,'r');

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tORF-fragmented reference FASTA : " + RNAseq_FASTA_file);

#............................................................................................................

# Setup array and counter for tracking fragment definition data.
fragments        = [];
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tProcessing restriction-digested genome file -> fragment coordinates.");

# Process digested FASTQ genome file, line by line.
for line in RNAseq_FASTA_data:
	first_char = line[:1];   # a '>' indicates a header line...
	if first_char == ">":
		# Line is header to FASTQ entry.
		line_parts             = string.split(string.strip(line));
		chrGenomeAndNum_string = line_parts[0];
		bp_coordinate_string   = line_parts[1];
		fragment_size_string   = line_parts[2];
		if len(line_parts) > 3:   # only usable lines should have a fourth space-delimited string in the header line.
			fragment_usable_string = line_parts[3];
			if fragment_usable_string[1] == "*":   # the second character of the final substring is '*' for useful fragments.
				# Fragment is usable, so the details should be placed into fragments structure.

				# split the chr string by '.' character, then trim off the first three characters ('chr') from the second substring.
				#   string has format of : ">Ca_a.chr1"
				genomeName_string,chrNum_string = chrGenomeAndNum_string.split(".");
				chr_num                         = int(float(chrNum_string.replace("chr","")));

				#   string has format of : "(9638..10115)"
				coordinates    = bp_coordinate_string.replace('(','').replace(')','').replace('..',' ').split();
				bp_start       = int(float(coordinates[0]));
				bp_end         = int(float(coordinates[1]));

				fragments.append([chr_num,bp_start,bp_end]);
RNAseq_FASTA_data.close();
#------------------------------------------------------------------------------------------------------------
# End of code section to parse restriction fragments from genome.
#============================================================================================================


print "###\t", time.clock() - t0, "seconds to parse restriction fragments from digested genome.";
t1 = time.clock();
print "### Starting read count data processing.";


#============================================================================================================
# Process genome repetitiveness file to determine repetitiveness max and average.
#------------------------------------------------------------------------------------------------------------

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tProcessing genome repetitiveness file -> max and average repetitiveness scores per fragment.");

# Look up chromosome name strings for genome in use.
#     Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/figure_definitions.txt"
figureDefinition_file  = workingDir + '/figure_definitions.txt';
figureDefinitionFile   = open(figureDefinition_file,'r');
figureDefinitionData   = figureDefinitionFile.readlines();

# Example lines in figureDefinition_file:
#     Chr  Use   Label   Name                         posX   posY   width   height
#     1    1     Chr1    Ca21chr1_C_albicans_SC5314   0.15   0.8    0.8     0.0625
#     2    1     Chr2    Ca21chr2_C_albicans_SC5314   0.15   0.7    *       0.0625
#     0    0     Mito    Ca19-mtDNA                   0.0    0.0    0.0     0.0

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
		if chr_num <> 0:
			chrName[int(float(chr_num))-1] = chr_name;
			with open(logName, "a") as myfile:
				myfile.write("\n\t\t\t\t\tChr" + str(chr_num) + " = " + chr_name);
figureDefinitionFile.close();

# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount;


#============================================================================================================
# Process repetitiveness score file into per chromosome lists.
#------------------------------------------------------------------------------------------------------------
# Open genome repetitiveness file.
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tOpen genome repetitiveness file : '" + datafile + "'");
data = open(datafile,'r');

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tGethering repetitiveness data for each chromosome.");
# Make per chromosome repetitive score lists.
repetData = [];
for x in chrName:
	repetData.append([]);
# Load data from file.
old_chr_name = '';
t0 = time.clock();
for line in data:
	# example lines from repetitiveness file:
	#     chromosome                    bpCoordinate    repetitivenessScore
	#     Ca21chr1_C_albicans_SC5314    2388924         123
	#     Ca21chr1_C_albicans_SC5314    2388925         135
	if line[0] == "#":
		print "### comment in repet file :'" + line.strip() + "'";
	else:
		line_parts = (line.strip()).split();
		chr_name   = line_parts[0];        # chr name of bp.
		position   = int(line_parts[1]);   # chr position of bp.
		repetScore = int(line_parts[2]);   # repetitiveness score at bp.

		if (chr_name <> old_chr_name):
			with open(logName, "a") as myfile:
				myfile.write("\n\t\t\t\t\tTime to process = " + str(time.clock() - t0) + " seconds.");
				myfile.write("\n\t\t\t\tLoading repetitiveness data for chr '" + str(chr_name) + "'");
			print '### Loading repetitiveness data for chr "' + str(chr_name) + '"';
			t0 = time.clock();

		# If chromosome string is being examined.
		if chr_name in chrName:
			# Identify the numerical placement of the chromosome in the chrName list.
			chrID_index = [i for i, x in enumerate(chrName) if x == chr_name];
			chr = chrID_index[0];

			# Append repetScore to per chr list.
			repetData[chr].append(repetScore);
		old_chr_name = chr_name;
data.close();


#============================================================================================================
# Output data per restriction fragment.
#------------------------------------------------------------------------------------------------------------
# fragments = [chr_num, bp_start, bp_end, data_sum, data_max, data_ave];
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tOutputting repetitiveness data of ORF-fragment genome.");
for current_fragment in fragments:
	chr_num         = current_fragment[0];
	bp_start        = current_fragment[1];
	bp_end          = current_fragment[2];
	chr_repet_data  = repetData[chr_num-1];
	repet_sum       = sum(chr_repet_data[(bp_start-1):bp_end]);
	repet_max       = max(chr_repet_data[(bp_start-1):bp_end]);
	fragment_length = bp_end - bp_start + 1;
	repet_ave       = repet_sum/float(fragment_length);
	# Output a line for each fragment: current_fragment = [chr_num, bp_start, bp_end]
	print str(chr_num) + '\t' + str(bp_start) + '\t' + str(bp_end) + '\t' + str(repet_max) + '\t' + str(repet_ave) + '\t' + str(fragment_length);
#------------------------------------------------------------------------------------------------------------
# End of code section to output information about fragments. 
#============================================================================================================


with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tWriting time to complete process to output files.")

print "### ", time.clock() - t0, "seconds to complete processing of pileup file and fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tTime to process = " + str(time.clock()-t0) )
	myfile.write("\n\t\t* 'scripts_genomes/genome_process_for_RNAseq.repetitiveness_2.py' completed. *\n")
