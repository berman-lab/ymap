# Input arguments: (Those with '[*]' at end are used here.)
#     2)  genome      : [String]: defines genome in use for project.             (ex. 'Ca_a') [*]
#     3)  workingDir  : [String]: Directory where links in system are stored.    (ex. '/home/bermanj/shared/links/) [*]
#
# Process input files:
#	1) Restriction-digested genome file.
#		*) Load usable fragment definitions into array : fragments[i][chr#,bpStart,bpEnd, data_count,data_max,data_ave]
#                                                                            [0   ,1      ,2    , 3         ,4       ,5       ]
#	2) Read counts for strain of interest dataset.
#		-) Find max read count on fragment.
#		-) Find average read count along fragment.
#
# Generate output file:
#	3) Output values for each fragment in a tab-delimited text file.
#		Each line contains information for one fragment = [chr_num,bp_start,bp_end, data_count,data_max,length]
#		0) chr_num  : Numerical chromosome identifier, defined for each genome in "figure_details.txt".
#		1) bp_start : Start bp coordinate along chromosome.
#		2) bp_end   : End bp coordinate along chromosome.
#		3) GC_ratio : Ratio of bases as 'G' or 'C' in fragment.
#	4) Comment lines in output begin with '###'.
#

import string, sys, re, time
workingDir  = sys.argv[1]
logName     = sys.argv[2]

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\n\t\t---------------------------------------------")
	myfile.write("\n\t\tLog of genome_process_for_RADseq.GC_bias_1.py")
	myfile.write("\n\t\t---------------------------------------------")

#============================================================================================================
# Process restriction-digested genome file.
#------------------------------------------------------------------------------------------------------------
# Example FASTQ header line.
#     >Ca_a.chr1 (9638..10115) (478bp) [*]
# FASTA entries with header lines ending in '[*]' are usable.

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tProcessing restriction-digested genome file.")

# Find name of genome FASTA file for species being examined.
#     Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/reference.txt"
reference_file = workingDir + '/reference.txt'
refFile        = open(reference_file,'r')
refFASTA       = refFile.read().strip()
refFile.close()

# Open restriction-digested genome FASTQ file.
ddRADseq_FASTA_file = workingDir + string.replace(refFASTA, '.fasta','.MfeI_MboI.fasta')
ddRADseq_FASTA_data = open(ddRADseq_FASTA_file,'r')

#............................................................................................................

# Setup array and counter for tracking fragment definition data.
fragments        = []
fragment_counter = 0

## Process digested FASTQ genome file, line by line.
while True:
	# Line pairs have the following structure.
	#    >Ca_a.chr1 (9638..10115) (478bp) [*]
	#    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
	line1 = ddRADseq_FASTA_data.readline()
	line2 = ddRADseq_FASTA_data.readline()
	if not line2:
		break  # EOF
	first_char = line1[:1];
	if first_char == ">":
		# Line is header to FASTQ entry.
		line_parts             = string.split(string.strip(line1))
		chrGenomeAndNum_string = line_parts[0]
		bp_coordinate_string   = line_parts[1]
		fragment_size_string   = line_parts[2]
		if len(line_parts) > 3:
			fragment_usable_string = line_parts[3]
			if fragment_usable_string[1] == "*":
				# Fragment is usable, so the details should be placed into fragments structure.
				#   chr number.
				#   fragment start bp.
				#   fragment end bp.

				# split the chr string by '.' character, then trim off the first three characters ('chr') from the second substring.
				#   string has format of : ">Ca_a.chr1"
				genomeName_string,chrNum_string = chrGenomeAndNum_string.split(".")
				chr_num                         = int(float(chrNum_string.replace("chr","")))

				#   string has format of : "(9638..10115)"
				coordinates = bp_coordinate_string.replace('(','').replace(')','').replace('..',' ').split()
				bp_start    = int(float(coordinates[0]))
				bp_end      = int(float(coordinates[1]))
				GC_ratio    = 0   # placeholder value.

				sequence = line2;
				G_count = sequence.count('G') + sequence.count('g')
				C_count = sequence.count('C') + sequence.count('c')
				T_count = sequence.count('T') + sequence.count('t')
				A_count = sequence.count('A') + sequence.count('a')
				if (float(G_count+C_count+T_count+A_count) == 0):
					GC_ratio = 0 
				else:
					GC_ratio = (G_count+C_count)/float(G_count+C_count+T_count+A_count)

				fragments.append([chr_num,bp_start,bp_end,GC_ratio])
				fragment_counter += 1
ddRADseq_FASTA_data.close()

# Put fragment counter into a general use variable.
numFragments = fragment_counter
#------------------------------------------------------------------------------------------------------------
# End of code section to parse restriction fragments from genome.
#============================================================================================================

print "### ", time.clock() - t0, "seconds to parse restriction fragments from digested genome."
t1 = time.clock()

print '### numFragments = ' + str(numFragments);
print '### Data from each fragment: [chrNum, bpStart, bpEnd, GC_ratio]'

#============================================================================================================
# Code section to output information about genome restriction fragments.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tOutputting GC-ratios of restriction-fragmented genome.")
for fragment in range(1,numFragments):
	# Output a line for each fragment.
	#     fragments[fragment-1] = [chr_num,bp_start,bp_end, GC_ratio]
	#     0) chr_num
	#     1) bp_start
	#     2) bp_end
	#     3) GC_ratio
	chr_num         = fragments[fragment-1][0]
	bp_start        = fragments[fragment-1][1]
	bp_end          = fragments[fragment-1][2]
	GC_ratio        = fragments[fragment-1][3]
	print str(chr_num) + '\t' + str(bp_start) + '\t' + str(bp_end) + '\t' + str(GC_ratio)

#------------------------------------------------------------------------------------------------------------
# End of code section to output information about fragments. 
#============================================================================================================

print "### ", time.clock() - t1, "seconds to output basic stats of each restriction fragment."
print "### ", time.clock() - t0, "seconds to complete processing of fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tTime to process = " + str(time.clock()-t0) )
	myfile.write("\n\t\t* 'py/genome_process_for_RADseq.GC_bias_1.py' completed. *")
