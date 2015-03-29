###
### Preprocesses child 'SNP_CNV_v1.txt' file that has been copied into the hapmap folder as
### "SNPdata_child.{entryNum}.txt".   Final format only has data lines for loci which are found
### in "SNPdata_hapmap.txt" file, and only data lines for fragments of the genome defined in the
### corresponding "haplotypeFragments.{entryNum}.txt" file.
### Output collumns:
###     0) chr_name      : Chromosome name string.
###     1) bp_coordinate : bp coordinate along chromosome.
###     2) baseCall_1    : majority allele base call.
###     3) baseCall_2    : secondary allele base call.
###

import string, sys, time;
genome     = sys.argv[1];
genomeUser = sys.argv[2];
project    = sys.argv[3];
user       = sys.argv[4];
hapmap     = sys.argv[5];
main_dir   = sys.argv[6];
entryNum   = sys.argv[7];
logName    = main_dir+"users/"+user+"/hapmaps/"+hapmap+"/process_log.txt";
inputFile1 = main_dir+"users/"+user+"/hapmaps/"+hapmap+"/SNPdata_parent.txt";
inputFile2 = main_dir+"users/"+user+"/hapmaps/"+hapmap+"/SNPdata_child."+entryNum+".txt";

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("*=======================================================================*\n");
	myfile.write("| Log of 'scripts_seqModules/scripts_hapmaps/hapmap.process_child.py'   |\n");
	myfile.write("*-----------------------------------------------------------------------*\n");
	myfile.write("| processing haplotype map entry "+entryNum+"\n");


#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+genomeUser+"/genomes/"+genome+"/";


#============================================================================================================
# Process genome LOH fragment file.
#------------------------------------------------------------------------------------------------------------
# Example fragment line.
#     >Ca_a.chr1 (9638..10115) (478bp) [*] a
with open(logName, "a") as myfile:
	myfile.write("|\tProcessing haplotype map entry '"+entryNum+"' genome LOH fragment file.\n");
# Open fragment file for haplotype map entry.
haplotypeMapEntry_file = main_dir+"users/"+user+"/hapmaps/"+hapmap+"/haplotypeFragments."+entryNum+".txt";
haplotypeMapEntry_data = open(haplotypeMapEntry_file,'r');
# Setup array and counter for tracking fragment definition data.
fragments        = [];
fragment_counter = 0;
## Process fragment file, line by line.
print '### Fragments found in haplotype entry:'
while True:
	# Line pairs have the following structure.
	#    >Ca_a.chr1 (9638..10115) (478bp) [*] a
	#    NULL
	line1 = haplotypeMapEntry_data.readline();
	line2 = haplotypeMapEntry_data.readline();
	if not line2:
		break  # EOF
	first_char = line1[:1];
	if first_char == ">":
		# Line is header to FASTA entry.
		line_parts             = string.split(string.strip(line1));
		chrGenomeAndNum_string = line_parts[0];
		bp_coordinate_string   = line_parts[1];
		fragment_size_string   = line_parts[2];
		homologID              = line_parts[4];
		# Fragment is usable, so the details should be placed into fragments structure.
		# split the chr string by '.' character, then trim off the first three characters ('chr') from the second substring.
		#   string has format of : ">Ca_a.chr1"
		genomeName_string,chrNum_string = chrGenomeAndNum_string.split(".");
		chr_num                         = int(float(chrNum_string.replace("chr","")));
		#   string has format of : "(9638..10115)"
		coordinates  = bp_coordinate_string.replace('(','').replace(')','').replace('..',' ').split();
		bp_start     = int(float(coordinates[0]));
		bp_end       = int(float(coordinates[1]));
		fragments.append([chr_num,bp_start,bp_end,homologID]);
		print '###\tfragment[' + str(fragment_counter) + '] = [' + str(chr_num) + ', ' + str(bp_start) + ', ' + str(bp_end) + ', ' + homologID + ']'
		fragment_counter += 1;
haplotypeMapEntry_data.close();
with open(logName, "a") as myfile:
	myfile.write("|\tThe genome LOH fragment definitions have been loaded.\n");
# Put fragment counter into a general use variable.
numFragments = fragment_counter;
print '### Number of fragments = '+str(numFragments);

#============================================================================================================
# Load FastaName from 'reference.txt' for genome in use.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tIdentifying name of reference FASTA file.\n");
reference_file = genomeDirectory + '/reference.txt';
refFile        = open(reference_file,'r');
FastaName      = refFile.read().strip();
refFile.close();
FastaName      = FastaName.replace(".fasta", "");


#============================================================================================================
# Look up chromosome name strings for genome in use; read in and parse genome 'figure_definitions.txt' file.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tLoad chromosome name strings from 'figure_definitions.txt' file.\n");
figureDefinition_file  = genomeDirectory + 'figure_definitions.txt';
figureDefinitionFile   = open(figureDefinition_file,'r');
figureDefinitionData   = figureDefinitionFile.readlines();
# Example lines in figureDefinition_file:
#     Chr  Use   Label   Name                         posX   posY   width   height
#     1    1     Chr1    Ca21chr1_C_albicans_SC5314   0.15   0.8    0.8     0.0625
#     2    1     Chr2    Ca21chr2_C_albicans_SC5314   0.15   0.7    *       0.0625
#     0    0     Mito    Ca19-mtDNA                   0.0    0.0    0.0     0.0


#============================================================================================================
# Determine the number of chromosomes of interest in genome.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tDetermining number of chromosomes of interest in genome.\n");
chrName_maxcount = 0;
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line));
	chr_num = line_parts[0];
	if chr_num.isdigit():
		chr_num    = int(float(line_parts[0]));
		chr_use    = int(float(line_parts[1]));
		chr_label  = line_parts[2];
		chr_name   = line_parts[3];
		if chr_num > chrName_maxcount:
			chrName_maxcount = chr_num;
figureDefinitionFile.close();


#============================================================================================================
# Gather name strings for chromosomes, in order.
#------------------------------------------------------------------------------------------------------------
# Pre-allocate chrName_array
chrName = [];
for x in range(0, chrName_maxcount+10):
	chrName.append([]);
with open(logName, "a") as myfile:
	myfile.write("|\tGathering name strings for chromosomes.\n");
figureDefinitionFile  = open(figureDefinition_file,'r');
chrCounter = 0;
chrNums    = [];
chrNames   = [];
chrLabels  = [];
chrShorts  = [];
for line in figureDefinitionData:
	line_parts = string.split(string.strip(line));
	chr_num = line_parts[0];
	if chr_num.isdigit():
		chr_num                        = int(float(line_parts[0]));
		chrNums.append(chr_num);
		chr_use                        = int(float(line_parts[1]));
		chr_label                      = line_parts[2];
		chrLabels.append(chr_label);
		chr_name                       = line_parts[3];
		chrNames.append(chr_name);
		chr_nameShort                  = chr_label;
		chrShorts.append(chr_nameShort);
		chrName[chrCounter] = chr_name;
		with open(logName, "a") as myfile:
			myfile.write("|\t\t" + str(chr_num) + " : " + chr_name + " = " + chr_nameShort + "\n");
		chrCounter += 1;
figureDefinitionFile.close();
# Put the chromosome count into a smaller name for later use.
chrCount = chrName_maxcount;
with open(logName, "a") as myfile:
	myfile.write("|\t\tMax chr string : " + str(chrCount) + "\n");


#============================================================================================================
# Process 'SNPdata_hapmap.txt' file to determine parent's initial SNP loci.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tProcessing datasetParent 'SNPdata_hapmap.txt' file => het loci.\n");
count            = 0;
old_chr          = 0;
fragment_found   = 0;
last_fragment    = 0;
current_fragment = 0;
log_count        = 0;
log_offset       = 0;

print '### InputFile1 = ' + inputFile1;
print '### Number of Chromosomes = ' + str(chrCount);
for x in range(0,chrCount):
	print '### \t' + str(x+1) + ' : ' + str(chrName[x]);
with open(logName, "a") as myfile:
	myfile.write("|\tProcessing parental SNP loci and determining which genome LOH fragment all loci belong to.\n");
# Process 'SNPdata_hapmap.txt' file, line by line.
data         = open(inputFile1,'r');
searchTarget = open(inputFile2,'r');
childLine    = chrName[0]+'\t1';
for line in data:
	# example lines from SNP pileup file:
	#   0) chr_name      : Chromosome name string.
	#   1) bp_coordinate : bp coordinate along chromosome.
	#   2) baseCall_1    : alphabetical first call of heterozgous locus.
	#   3) baseCall_2    : alphabetical second call of heterozygous locus.
	if line[0] != "#":
		count += 1;
		parentLine = string.strip(line);
		parentLine = parentLine.split('\t');
		P_chr_name = parentLine[0];        # chr name of bp.			: Ca21chrR_C_albicans_SC5314
		P_position = int(parentLine[1]);   # chr position of bp.		: 2286371
		P_allele1  = parentLine[2];        # alphabetical 1' heterozygous allele.
		P_allele2  = parentLine[3];        # alphabetical 2' heterozygous allele.
		#=====================================
		# Find fragment bin which holds data.
		#-------------------------------------
		# Identify which chromosome this data point corresponds to.
		chr = 0;
		for x in range(0,chrCount):
			if (chrNums[x] != 0):
				if chrName[x] == P_chr_name:
					chr = x+1;
		# Display status updates to log file.
		if old_chr != chr:
			print '### chr change : ' + str(old_chr) + ' -> ' + str(chr);
			with open(logName, "a") as myfile:
				myfile.write("|\t    " + str(old_chr) + " -> " + str(chr) + " = " + P_chr_name + "\n");
				myfile.write("|\t\t1........01........01........01........01........01........01........01........01........01........0\n");
		if chr!=0:
			# Reset for each new chromosome examined.
			if old_chr != chr:
				if log_offset != 0:
					log_offset_string = " "*((log_offset)%100);
					with open(logName, "a") as myfile:
						myfile.write("|\t\t" + log_offset_string);
				count            = 1;
				fragment_found   = 0;
				current_fragment = 0;
			# If (fragment_found == 0), look to see if current coordinate matches any defined fragments.
			if fragment_found == 0:
				for frag in range(1,numFragments+1):
					# Check if current coordinate is consistent with this fragment : fragments[frag-1] = [chr_num, bp_start, bp_end, homologID]
					if chr == fragments[frag-1][0] and P_position >= fragments[frag-1][1] and P_position <= fragments[frag-1][2]:
						fragment_found    = 1;
						current_fragment  = frag;
						break;
			# If (fragment_found == 1), add current bp coordinate data to fragment data.
			if fragment_found == 1:
				# display status updates to log file.
				if (count%50000 == 0):
					log_count  += 1;
					log_offset += 1;
					if (log_count%100) == 0:
						if (log_count-1) == 0:
							with open(logName, "a") as myfile:
								myfile.write("|\t\t");
						else:
							with open(logName, "a") as myfile:
								myfile.write(" " + str(log_count-1) + "\n");
					with open(logName, "a") as myfile:
						myfile.write(".");
				# If current coordinate is at/after end of a fragment, update fragment_found to 0.
				if P_position >= fragments[current_fragment-1][2]:
					fragment_found = 0;
			#===============================================================================================================
			# Get line corrosponding to parental het SNP in child 'SNPdata_child.{entryNum}.txt' file.
			#--------------------------------------------------------------------------------------------------------------- 
			searchString = P_chr_name + "\t" + str(P_position);
			result = "";
			## Barely functional, very slow slow equivalent of GREP, returning only one result line.
			##    Very slow because it restarts search position for each search.
			##    This can be improved dramatically if the search target file is already sorted.
			#for lineChild in open(inputFile2,'r'):
			#	if searchString in lineChild:
			#		result = lineChild
			#		break
			# Fast line search, relies upon presorted search target file.
			if len(childLine) == 0:
				childLine       = searchTarget.readline();
			childLine_parts     = childLine.split('\t');
			C_chr_name          = childLine_parts[0];
			C_position          = int(childLine_parts[1]);
			C_chr = 0;
			for x in range(0,chrCount):
				if (chrNums[x] != 0):
					if chrName[x] == C_chr_name:
						C_chr = x+1;
			P_chr = chr;
			while P_chr > C_chr:	# WORKING: this section jumps through the child lines of chromosomes with no parent lines.
				childLine       = searchTarget.readline();
				childLine       = string.strip(childLine);
				childLine_parts = childLine.split('\t');
				C_chr_name      = childLine_parts[0];
				C_position      = int(childLine_parts[1]);
				for x in range(0,chrCount):
					if (chrNums[x] != 0):
						if chrName[x] == C_chr_name:
							C_chr = x+1;
			while P_chr == C_chr and P_position > C_position:   # WORKING: this section jumps through the child lines until the correct chromosome and coordinate is reached.
				childLine       = searchTarget.readline();
				childLine       = string.strip(childLine);
				childLine_parts = childLine.split('\t');
				C_chr_name      = childLine_parts[0];
				C_position      = int(childLine_parts[1]);
				for x in range(0,chrCount):
					if (chrNums[x] != 0):
						if chrName[x] == C_chr_name:
							C_chr = x+1;
			if P_chr > C_chr:
				result          = "";
			elif P_chr == C_chr and P_position > C_position:
				result          = "";
			elif P_chr == C_chr and P_position < C_position:
				result          = "";
			elif P_chr == C_chr and P_position == C_position:
				result          = childLine;

			if result == "":
				locus_phase = 10;  # no phase information, due to no data at coodinate in child dataset.
			else:
				C_chr_name    = childLine_parts[0];   # chr_name_string
				C_position    = childLine_parts[1];   # bp_coordinate
				C_countTot    = childLine_parts[2];   # total_reads
				C_refBase     = childLine_parts[3];   # reference_base
				C_countA      = childLine_parts[4];   # A reads
				C_countT      = childLine_parts[5];   # T reads
				C_countG      = childLine_parts[6];   # G reads
				C_countC      = childLine_parts[7];   # C reads
				# Child locus should be homozygous, or at least unbalanced, so figure out which allele is most common.
				alleleData    = [(int(float(C_countA)), 'A'), (int(float(C_countT)), 'T'), (int(float(C_countG)), 'G'), (int(float(C_countC)), 'C')];
				sortedAlleles = sorted(alleleData, key=lambda alleleDatum: alleleDatum[0]); # sort alleles by copy number.
				C_allele      = sortedAlleles[3][1];
				# Retrieve the expected allele from the genome LOH fragments in the haplotype map.
				#     fragments[current_fragment] = [chr_num,bp_start,bp_end,homologID]
				homologID     = fragments[current_fragment-1][3];
				# Determine phasing of SNP locus by comparing found child allele to haplotype map entry expectation.
				#     no phase info => 0
				#     correct phase => 1
				#     wrong phase   => 2 
				if (fragment_found == 0): # this coordinate was not associated with a LOH fragment.
					locus_phase = 11;          # no phase information, due to coordinate not matching a LOH fragment in child dataset.
				else:
					if (C_allele == P_allele1):
						if (homologID == "a") or (homologID == "A"):
							locus_phase = 0;   # correct phase.
						else:
							locus_phase = 1;   # wrong phase.
					elif (C_allele == P_allele2):
						if (homologID == "a") or (homologID == "A"):
							locus_phase = 1;   # wrong phase.
						else:
							locus_phase = 0;   # correct phase.
					else:
						locus_phase = 12;      # no phase information, due to surprise allele in child dataset.
			# Output 'SNPdata_hapmap.txt' line with new column containing most recent haplotype entry phase.
			print string.strip(line) + '\t' + str(locus_phase) #+ '\t(' + str(current_fragment) + ';' + homologID + ')';
		# Reset old_chr/last_fragment to current coordinate chromosome before moving to next line in pileup. 
		old_chr       = chr;
		last_fragment = current_fragment;
data.close();
searchTarget.close();


#============================================================================================================
# Final log file output.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("|\tTime to process = " + str(time.clock()-t0) + "\n");
	myfile.write("*-----------------------------------------------------------------------*\n");
	myfile.write("| End of 'scripts_seqModules/scripts_hapmaps/hapmap.process_child.py'   |\n");
	myfile.write("*=======================================================================*\n");
