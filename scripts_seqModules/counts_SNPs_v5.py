# Processing SAM/BAM files for putative SNPs.
#------------------------------------------------------------------------------------------------------------
# Generate pileup file using samtools:
#	"samtools pileup -f genome.fasta my_file.bam | awk '{print $1 " " $2 " " $3 " " $4 " " $5}' > output_file.pileup"
# Generate putative SNPs list using this script:
# 	"python counts_SNPs_v3.py FH1.pileup > FH1_putative_SNPs_v#.txt"
#============================================================================================================

import string, sys, re

# python 2 : my_file = file(sys.argv[1],'r').xreadlines()
# python 3
my_file = open(sys.argv[1]);


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
				if astr[i] == "A" or astr[i] == "T" or astr[i] == "G" or astr[i] == "C" or astr[i] == "." or astr[i] == ",": result += astr[i]
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
	line                              = i.strip().split();
	chrom                             = line[0]				# chromosome label for locus.
	pos                               = line[1]				# coordinate for locus (in bp).
	ref_base                          = line[2]				# reference base at this locus.
	total                             = line[3]				# total count of reads at locus.
	if (len(line) > 4):
		reads                     = line[4].upper()			# string defining locus => capitalized.
		reads_noStartEnd          = dump_startend(reads)		# locus string without indels.
		reads_noIndels_noStartEnd = dump_indels(reads_noStartEnd)	# locus string without indels or end/start/quality.
		A                         = len(re.findall("A", reads_noIndels_noStartEnd))
		T                         = len(re.findall("T", reads_noIndels_noStartEnd))
		G                         = len(re.findall("G", reads_noIndels_noStartEnd))
		C                         = len(re.findall("C", reads_noIndels_noStartEnd))
		ref_count                 = len(re.findall("\.", reads_noIndels_noStartEnd)) + len(re.findall("\,", reads_noIndels_noStartEnd))
	else:
		A                         = 0;
		T                         = 0;
		G                         = 0;
		C                         = 0;
		ref_count                 = 0;
		# count of reads identical to reference at this locus.

	#print chrom + '\t' + pos + '\t' + ref_base + '\t' + str(A) + '\t' + str(T) + '\t' +  str(G) + '\t' +  str(C) + '\t' +  str(ref_count)
	# Adds reference base count to appropriate counter.
	# Without this, the apparent reads would only account for variations from reference, not the
	#    reference itself.   "...,...,,,...T.." would be interpreted as a single base (T) seen, instead of
	#    two (T and ref).
	# if samtools output does not contain '.,' characters for matching to reference, then ref_base = 'N'
	#    and no correction is needed for previous 'ATCG' counts.
	if ref_base == "A":
		A = ref_count
	elif ref_base == "T":
		T = ref_count
	elif ref_base == "C":
		C = ref_count
	elif ref_base == "G":
		G = ref_count

	# boolean interpretation of alternate bases from reference present in reads for locus.
	# isA+isT+isG+isC > 1 when more than one base is seen at this locus.
	isA = isT = isG = isC = 0
	if A != 0: isA = 1
	if T != 0: isT = 1
	if G != 0: isG = 1
	if C != 0: isC = 1

	if isA+isT+isG+isC > 1:		# Only deal with loci where more than one base is seen.
		print(chrom + '\t' + pos + '\t' + ref_base + '\t' + str(A) + '\t' + str(T) + '\t' +  str(G) + '\t' +  str(C))

