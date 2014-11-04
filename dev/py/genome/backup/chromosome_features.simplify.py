# py/genome/chromosome_features.simplify.py 
#
# Input Arguments
#   1) user       : 'darren'
#   2) genome     : 'Candida_albicans_SC5314_etc'
#   3) main_dir   : '/home/bermanj/shared/links/'
#   4) logName    :
#
# Generate output file:
#   1) a simplified pileup file containing average read counts per RADseq feature.   [chr_num,bp_start,bp_end, data_ave]
#       0) chr_num   : Numerical chromosome identifier, defined for each genome in "figure_details.txt".
#       1) bp_start  : Start bp coordinate along chromosome.
#       2) bp_end    : End bp coordinate along chromosome.
#

import string, sys, re, time
user       = sys.argv[1]
genome     = sys.argv[2]
main_dir   = sys.argv[3]
logName    = sys.argv[4]

t0 = time.clock()

with open(logName, "a") as myfile:
	myfile.write("\t\t\t*================================================================*\n")
	myfile.write("\t\t\t| Log of 'chromosome_features.simplify.py'                       |\n")
	myfile.write("\t\t\t*----------------------------------------------------------------*\n")

#============================================================================================================
# Find location of genome being used.
#------------------------------------------------------------------------------------------------------------
genomeDirectory = main_dir+"users/"+user+"/genomes/"+genome+"/"


#============================================================================================================
# Load file description information from 'expression.txt' for genome in use.
#------------------------------------------------------------------------------------------------------------
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tLoading chromosome_features file information from 'expression.txt'.\n")
expression_file = genomeDirectory + 'expression.txt'
expressionFile  = open(expression_file,'r')
exp_header_rows = int(float(expressionFile.readline().strip()))
exp_chrId_col   = int(float(expressionFile.readline().strip()))
exp_startBP_col = int(float(expressionFile.readline().strip()))
exp_endBP_col   = int(float(expressionFile.readline().strip()))
expressionFile.close()


#============================================================================================================
# Process chromosome feature file.
#------------------------------------------------------------------------------------------------------------
## print "### Processing feature file."
# header lines from "chromosome_features.txt" defined by [exp_header_rows].
# example line from "chromosome_features.txt", defined by [exp_chrId_col, exp_startBP_col, exp_endBP_col].
feature_file     = genomeDirectory + "chromosome_features.txt"
feature_dataFile = open(feature_file,'r')
# Setup array and counter for tracking feature definition data.
features        = []
feature_counter = 0
# Skip through header lines.
for x in range(0,exp_header_rows):
	discard_header_line = feature_dataFile.readline()
# Process remaining lines in 'chromosome_features.txt' file, line by line.
for dataLine in feature_dataFile:
	line_parts   = dataLine.split("\t")
	if len(line_parts) > max([exp_chrId_col-1,exp_startBP_col-1,exp_endBP_col-1]):
		if len(line_parts[exp_chrId_col-1]) > 0:
			chrID_string = line_parts[exp_chrId_col-1]
			startBP_int  = int(float(line_parts[exp_startBP_col-1]))
			endBP_int    = int(float(line_parts[exp_endBP_col-1]))
			# Example line
			#	chrID_string : 'Ca21chr3_C_albicans_SC5314'
			#	startBP_int  : 394543
			#	endBP_int    : 390791

			# if startBP_int and endBP_int are not in the proper order, swap their values.
			if endBP_int < startBP_int:
				temp_int    = startBP_int
				startBP_int = endBP_int
				endBP_int   = temp_int

			# Add feature entry to list.
			features.append([chrID_string,startBP_int,endBP_int])
			feature_counter += 1
feature_dataFile.close()
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tThe chromosome feature file for the genome has been loaded.\n")
# Put feature counter into a general use variable.
numFeatures = feature_counter

# Sort collected data by start coordinate, then by chromosome name.  Because sorts are 'guaranteed to be stable', this will result in the proper overall arrangement.
## print "### Sort 1 of feature list."
features = sorted(features, key = lambda x: x[1] )
## print "### Sort 2 of feature list."
features = sorted(features, key = lambda x: x[0] )
#------------------------------------------------------------------------------------------------------------
# End of code section to parse chromosome features file.
#============================================================================================================


#============================================================================================================
# Code section to output simplified and sorted chromosome feature information.
#------------------------------------------------------------------------------------------------------------
## print "### Output of new feature list."
with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tOutput simplified and sorted feature list.\n")
for feature in range(1,numFeatures):
	# Output a line for each feature.
	#     features[feature-1] = [chr_num,bp_start,bp_end, aveDepth]
	#     0) chromosome_name
	#     1) bp_start
	#     2) bp_end
	chrID_string = features[feature-1][0]
	startBP_int  = features[feature-1][1]
	endBP_int    = features[feature-1][2]
	print chrID_string + '\t' + str(startBP_int) + '\t' + str(endBP_int)
#------------------------------------------------------------------------------------------------------------
# End of code section to output information about features. 
#============================================================================================================


with open(logName, "a") as myfile:
	myfile.write("\t\t\t|\tNew simplified and sorted feature file has been generated.\n")
	myfile.write("\t\t\t*----------------------------------------------------------------*\n")
	myfile.write("\t\t\t| End of Log from 'chromosome_features.simplify.py '             |\n")
	myfile.write("\t\t\t*================================================================*\n")
