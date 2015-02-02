# Input arguments:
#	1) user    : Name of user account installing the genome file.
#	2) genome  : Name of the genome.
#	3) logfile : Path and file name of output log file.
#
# Output:
#

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
# Process used chromosomes into k-mer array.
#------------------------------------------------------------------------------------------------------------

# Open reformatted FASTA file.
FASTA_file = workingDir + refFASTA;
FASTA_data = open(FASTA_file,'r');

# Process reformatted FASTA file, entry by entry.
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
		line_parts             = string.split(string.strip(line1))
		chromosome_name        = line_parts[0]
		bp_coordinate_string   = line_parts[1]
		fragment_size_string   = line_parts[2]
		if len(line_parts) > 3:



#============================================================================================================
# 
#------------------------------------------------------------------------------------------------------------


#============================================================================================================
# 
#------------------------------------------------------------------------------------------------------------


#............................................................................................................

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tOpen genome repetitiveness file.")

# Open genome repetitiveness file.
# datafile    = workingDir + genome + '_repetitiveness.txt'
print '### InputFile = ' + inputFile
#datafile      = workingDir + inputFile
datafile      = inputFile;
data          = open(datafile,'r')

#............................................................................................................

count            = 0
old_chr          = 0
fragment_found   = 0
last_fragment    = 0
current_fragment = 0
log_offset       = 0

print '### Number of Chromosomes = ' + str(chrCount)
for x in range(0,chrCount):
	print '### \t' + str(x) + ') ' + str(chrName[x])

print "###" + str(numFragments)

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tGethering repetitiveness data for each fragment.")

# Process repetitiveness file, line by line.
for line in data:
	# example lines from SNP pileup file:
	#     chromosome                    bpCoordinate    repetitivenessScore
	#     Ca21chr1_C_albicans_SC5314    2388924         123
	#     Ca21chr1_C_albicans_SC5314    2388925         135

	if line[1] == "#":
		print '### comment in repet file : "' + string.strip(line) + '"'
	else:
		count += 1
		line_parts = string.split(string.strip(line))
		chr_name   = line_parts[0]   # chr name of bp.
		position   = line_parts[1]   # chr position of bp.
		repetScore = line_parts[2]   # repetitiveness score at bp.

		# Attempt to match up current data line with pre-determined restriction fragments.
		found = 0

		# Identify which chromosome this data point corresponds to.
		chr = 0
		for x in range(0,chrCount):
			# print str(chrName[x])
			if chrName[x] == chr_name:
				chr = x+1

		# Convert to integers.
		pos_point  = int(position)   
		data_point = float(repetScore)

		if old_chr != chr:
			print '### chr change : ' + str(old_chr) + ' -> ' + str(chr)
			with open(logName, "a") as myfile:
				myfile.write("\n\t\t\t\t    chr : " + str(old_chr) + " -> " + str(chr))
				myfile.write("\n\t\t\t\t\t1........01........01........01........01........01........01........01........01........01........0")

		if chr!=0:
			# Reset for each new chromosome examined.
			if old_chr != chr:
				if log_offset != 0:
					log_offset_string = " "*((log_offset)%100)
					with open(logName, "a") as myfile:
						myfile.write("\n\t\t\t\t\t" + log_offset_string)
				count = 1
				fragment_found = 0
				current_fragment = 0

			# If (fragment_found == 0), look to see if current coordinate matches any defined fragments.
			if fragment_found == 0:
				for frag in range(current_fragment,numFragments):
					# Check if current coordinte is consistent with this fragment : fragments[frag-1] = [chr_num,bp_start,bp_end,
					#       data_count,data_max,ave_read_count, repetDataCount,repetMax,repetAve]
					if chr == fragments[frag-1][0] and pos_point >= fragments[frag-1][1] and pos_point <= fragments[frag-1][2]:
						fragment_found   = 1
						current_fragment = frag
						break

			# If (fragment_found == 1), add current bp coordinate data to fragment data.
			if fragment_found == 1:
				if last_fragment != current_fragment:
					log_offset += 1
					if ((current_fragment-1)%100) == 0:
						if (current_fragment-1) == 0:
							with open(logName, "a") as myfile:
								myfile.write("\n\t\t\t\t\t")
						else:
							with open(logName, "a") as myfile:
								myfile.write(" " + str(current_fragment-1) + "\n\t\t\t\t\t")
					with open(logName, "a") as myfile:
						myfile.write(".")

				# Adds current coordinate repetitiveness score to fragment total count : fragments[frag-1] = [chr_num,bp_start,bp_end,
				#       repetDataCount,repetMax,repetAve]
				fragments[current_fragment-1][3] += data_point

				# If current coordinate read count is highest so far for fragment, update max : fragments[frag-1] = [chr_num,bp_start,bp_end,
				#       repetDataSum,repetMax,repetAve]
				if data_point > fragments[current_fragment-1][4]:
					fragments[current_fragment-1][4] = data_point

				# If current coordinate is at (or after) the end of a fragment, update fragment_found to 0 : fragments[frag-1] = [chr_num,bp_start,bp_end,
				#       repetDataSum,repetMax,repetAve]
				if pos_point >= fragments[current_fragment-1][2]:
					fragment_found = 0

		# Reset old_chr to current coordinate chromosome before moving to next line in pileup. 
		old_chr = chr
		last_fragment = current_fragment

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tCalculating average read coverage per fragment.")

# Calculate average read coverage per each fragment.
for fragment in range(1,numFragments):
	# Calculate average read coverage for this fragment.
	#       fragments[fragment-1] = [chr_num,bp_start,bp_end, repetDataSum,repetMax,repetAve]
	fragments[fragment-1][5] = fragments[fragment-1][3]/float(fragments[fragment-1][2]-fragments[fragment-1][1]+1)
#------------------------------------------------------------------------------------------------------------
# End of code section to parse repetitiveness max/average data.
#============================================================================================================


print "### ", time.clock() - t1, "seconds to parse genome repetitiveness file."
t2 = time.clock()

print '### Number of fragments = ' + str(numFragments)
print '### Data from each fragment: [chrNum, bpStart, bpEnd, Max, Ave, Length]'

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

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t\tWriting time to complete process to output files.")

print "### ", time.clock() - t2, "seconds to output basic stats of each restriction fragment."
print "### ", time.clock() - t0, "seconds to complete processing of pileup file and fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tTime to process = " + str(time.clock()-t0) )
	myfile.write("\n\t\t* 'py/genome_process_for_standard_bins.repetitiveness_1.py' completed. *\n")
