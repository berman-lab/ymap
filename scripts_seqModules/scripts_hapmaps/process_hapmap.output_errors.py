#
# usage:
#	python process_hapmap.output_errors.py [raw_hapmap_file] > [errors_only_file]
#
# The [raw_hapmap_file] is probably "SNPdata_parent.txt"
# The [errors_only_file] is up to the user.
#

def process_HapmapLine(entry_line):
	# Process 'putative_SNPs_v4.txt' file line.
	# example lines:
	#       chromosome                   coord   ref   A    T     G   C
	#       Ca21chrR_C_albicans_SC5314	1748310	A	G	1
	#       Ca21chrR_C_albicans_SC5314	1753108	C	G	1
	#       Ca21chrR_C_albicans_SC5314	1759210	G	A	1
	#       Ca21chrR_C_albicans_SC5314	1772100	G	A	1
	parent_line = string.strip(entry_line)
	parent_line = parent_line.split('\t')
	P_chr_name  = parent_line[0]   # chr name of bp.          : Ca21chrR_C_albicans_SC5314
	P_position  = parent_line[1]   # chr position of bp.      : 2286371
	P_homologA  = parent_line[2]   # 
	P_homologB  = parent_line[3]   # 
	P_status    = parent_line[4]   # 0 = good; 1 = need to flip; 10/11/12 = bad.
	return P_chr_name,P_position,P_homologA,P_homologB,P_status

import string, sys, time

inputFileName     = sys.argv[ 1]
data_P = open(inputFileName,"r")
line_P = data_P.readline()
error_endOfFile = False
while (error_endOfFile == False):
	first_char = line_P[:1]
	if first_char <> "#":
		P_chrName,P_position,P_homologA,P_homologB,P_status = process_HapmapLine(line_P)
		if (int(P_status) > 1):
			print P_chrName+"\t"+P_position+"\t"+P_homologA+"\t"+P_homologB+"\t"+P_status

	line_P = data_P.readline()
	if not line_P: # EOF 1
		error_endOfFile = True
		break
data_P.close()
