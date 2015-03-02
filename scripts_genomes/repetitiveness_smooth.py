# Input arguments: (Those with '[*]' at end are used here.)
# Process genome repetitiveness file by restriction fragment definitions.
#	1) test.repetitiveness.txt		(Ymap_root/users/[user]/genomes/[genome]/)
#	2) test.MfeI_MboI.fasta			(Ymap_root/users/[user]/genomes/[genome]/)
# Output fragment repetitiveness data.
#	1) test.repetitiveness.MfeI_MboI.txt	(Ymap_root/users/[user]/genomes/[genome]/)
#
import string, sys, re, time, math, numpy;
userName        = sys.argv[1];
genomeName      = sys.argv[2];
main_dir        = sys.argv[3];
logName         = sys.argv[4];
Smoothing_sigma = int(sys.argv[5]);

#------------------------------------------------------------------------------------------------------------
def gauss(n,sigma):
	r = range(-int(n/2),int(n/2)+1);
	gauss_list    = [1 / (sigma * math.sqrt(2*math.pi)) * math.exp(-float(x)**2/(2*sigma**2)) for x in r];
	gauss_list[:] = [x*sum(gauss_list) for x in gauss_list];
	return gauss_list; # ensure gaussian kernal totals to one.
#------------------------------------------------------------------------------------------------------------

t0 = time.clock();
with open(logName, "a") as myfile:
	myfile.write("\t\t\t*================================================================*\n");
	myfile.write("\t\t\t| Log of 'scripts_genomes/repetitiveness_smooth.py'              |\n");
	myfile.write("\t\t\t*----------------------------------------------------------------*\n");


#============================================================================================================
# Load genome definition files to determine which chromosomes are to be examined.
#------------------------------------------------------------------------------------------------------------
# Look up chromosome name strings for genome in use.
#     Read in and parse : "links_dir/main_script_dir/genome_specific/[genome]/figure_definitions.txt"
workingDir             = main_dir + 'users/' + userName + '/genomes/' + genomeName + '/';
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
datafile      = workingDir + genomeName + ".repetitiveness.txt";
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
	if (line[0] == "#"):
		print "### comment in repet file :'" + line.strip() + "'";
	elif (line.strip() == ''):
		print "###";
	else:
		line_parts = (line.strip()).split();
		chr_name   = line_parts[0];          # chr name of bp.
		position   = int(line_parts[1]);     # chr position of bp.
		repetScore = float(line_parts[2]);   # repetitiveness score at bp.

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
# Output smoothed chromosome repetitiveness data.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\n\t\t\tGaussian smoothing repetitiveness data.");
	myfile.write("\n\t\t\tGaussian smoothing kernel sigma = " + str(int(Smoothing_sigma)));
	myfile.write("\n\t\t\tGaussian smoothing kernel width = " + str(int(Smoothing_sigma*3+1)));
# Generate Gaussian kernel filter: gauss(width,sigma);
kernel_sigma = Smoothing_sigma;
kernel_width = kernel_sigma*3+1;
filter = gauss(kernel_width,kernel_sigma);
print "### Gaussian smoothed repetitiveness being calculated.";
for chr in range(0, chrName_maxcount):
	if len(chrName[chr]) > 0:
		print "### Smoothing " + str(chrName[chr]) + " for output.";
		# Convolve repetitiveness data.
		smoothed_chr = numpy.convolve(repetData[chr],filter).tolist();
		# Trim ends of length [(kernel_width-1)/2] from smoothed_chr.
		for i in range((kernel_width-1)/2):
			del smoothed_chr[0];
			del smoothed_chr[-1];
		# Output smoothed repetitiveness data.
		#   example lines from repetitiveness file:
		#     chromosome                    bpCoordinate    repetitivenessScore
		#     Ca21chr1_C_albicans_SC5314    2388924         123
		#     Ca21chr1_C_albicans_SC5314    2388925         135
		for pos in range(len(smoothed_chr)):
			print str(chrName[chr]) + "\t" + str(pos+1) + "\t" + str(smoothed_chr[pos]);
#------------------------------------------------------------------------------------------------------------
# End of code section to output Gaussian smoothed repetitiveness data. 
#============================================================================================================


print "### ", time.clock() - t0, "seconds to complete processing of pileup file and fragment definitions."

with open(logName, "a") as myfile:
	myfile.write("\n\t\t\t| Time to process = " + str(time.clock()-t0) + "\n")
	myfile.write(  "\t\t\t*----------------------------------------------------------------*\n");
	myfile.write(  "\t\t\t| 'scripts_genomes/repetitiveness_smooth.py' completed.          |\n");
	myfile.write(  "\t\t\t*================================================================*\n");
