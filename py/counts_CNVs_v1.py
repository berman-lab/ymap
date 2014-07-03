# Processing SAM/BAM files for CNV information.
#------------------------------------------------------------------------------------------------------------
# Generate pileup file using samtools:
#	"samtools pileup my_file.bam -s | awk '{print $1 " " $2 " " $3 " " $4 " " $5 " " $7}' > output_file.pileup"
# Generate putative SNPs list using this script:
# 	"python counts_CNVs_v3.py FH1.pileup > FH1_putative_CNVs_v#.txt"
#============================================================================================================

import string, sys, re
my_file = file(sys.argv[1],'r').xreadlines()

#------------------------------------------------------------------------------------------------------------
# find_maxQuality(astr) determines the maximum quality from a quality string, converts it to a numeric value.
def find_maxQuality(astr):
	maxQuality = 0
	for i in range(0, len(astr)):
		curQuality = ord(astr[i])
		if curQuality > maxQuality:
			maxQuality = curQuality
	return maxQuality
#------------------------------------------------------------------------------------------------------------
# find_aveQuality(astr) determines the maximum quality from a quality string, converts it to a numeric value.
def find_aveQuality(astr):
	totQuality = 0
	totCount   = 0;
	for i in range(0, len(astr)):
		curQuality = ord(astr[i])
		totQuality = totQuality+curQuality
	if totCount == 0:
		result = 0;
	else:
		result = totQuality/totCount
	return result

#------------------------------------------------------------------------------------------------------------
for i in my_file:	# process pileup file line by line.
	line     = string.split(string.strip(i), ' ')
	chrom    = line[0]			# chromosome label for locus.
	pos      = line[1]			# coordinate for locus (in bp).
	ref_base = line[2]			# reference base at this locus.
	total    = line[3]			# total count of reads at locus.
	quality  = line[4]			# 'max mapping quality'

	print chrom + '\t' + pos + '\t' + total
