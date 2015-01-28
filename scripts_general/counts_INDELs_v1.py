# Processing SAM/BAM files for putative SNPs.
#------------------------------------------------------------------------------------------------------------
# Generate pileup file using samtools:
#	"samtools pileup -f CaSC5314_v21.fasta my_file.bam -i > my_indels.pileup"
#	(chromosome name; coordinate; base; read count; reads; read quality)"
# Generate putative SNPs list using this script:
# 	"python counts_INDELs_v1.py my_indels.pileup > FH1_putative_INDELs_v#.txt"
#============================================================================================================

import string, sys, re
my_file = file(sys.argv[1],'r').xreadlines()

#------------------------------------------------------------------------------------------------------------
# dump_indels(astr) removes any insertions or deletions from the read base data in the pileup.
#	'\+[0-9]+[ATCGNatcgn]+' : indicates an insertion.
#	'\-[0-9]+[ATCGNatcgn]+' : indicates a deletion.
def dump_indels(astr):
        result = ""
        blackout = 0
        for i in range(0, len(astr)):
                if astr[i] == "+" or astr[i] == "-":
                        start = i+1
                        val = ""
                        while astr[start] > '0' and astr[start] <= '9':
                                val += astr[start]
                                start += 1
                        if val == "":
                                blackout = 1
                        else:
                                blackout = int(val) + 1
                else:
                        if blackout != 0:
                                blackout -= 1
                        else:
                                if astr[i] == "A" or astr[i] == "T" or astr[i] == "G" or astr[i] == "C": result += astr[i]
        return result

#------------------------------------------------------------------------------------------------------------
# dump_startend(astr) removes any marks indicating a start or end of a read segment from the read base data.
#	'\$'  : indicates the start of a read.
#	'\^.' : indicates the end of a read, with a single character describing quality of that read.
def dump_startend(astr):
	result   = ""
	lastchar = ""
	oldchar  = ""
	for i in range(0, len(astr)):
		newchar = astr[i]
		if newchar == "$":
			oldchar = ""
		elif oldchar == "^":
			oldchar = ""
		elif newchar == "^":
			oldchar = "^"
		else:
			result += newchar
			oldchar = newchar
	return result        

#------------------------------------------------------------------------------------------------------------
for i in my_file:	# process pileup file line by line.
	line                      = string.split(string.strip(i), ' ')
	chrom                     = line[0]				# chromosome label for locus.
	pos                       = line[1]				# coordinate for locus (in bp).
	total                     = line[3]				# total count of reads at locus.
	reads                     = string.upper(line[4])		# string defining locus => capitalized.
	reads_noStartEnd          = dump_startend(reads)		# locus string without indels.
	inserts                   = len(re.findall("\+", reads_noStartEnd))
	deletions                 = len(re.findall("\-", reads_noStartEnd))

	if inserts+deletions > 0:		# Only deal with loci with an INDEL.
		print chrom + '\t' + pos + '\t' + str(total) + '\t' + str(inserts) + '\t' + str(deletions)
