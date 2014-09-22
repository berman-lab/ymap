# Input Arguments
#   1) user     : 'darren'
#   2) hapmap   : 'test_Ca'
#   4) main_dir : '/home/bermanj/shared/links/'
#	5) logName  :
#
# Example input line:
#     2(b[1:495225] a[501338:3188548]); 3(a[1:1315157] b[1335699:2232035]);
#
# Example outputlines:
#     >Ca_haploid.b.chr2 (1..495225) (1091bp) [*]
#     NULL
#     >Ca_haploid.a.chr2 (501338..3188548) (2687211bp) [*]
#     NULL
#     >Ca_haploid.a.chr3 (1..1315157) (1315157bp) [*]
#     NULL
#     >Ca_haploid.b.chr3 (1335699..2232035) (896337bp) [*]

import string, sys
user       = sys.argv[1]
hapmap     = sys.argv[2]
main_dir   = sys.argv[3]
hapmapDir  = main_dir+"users/"+user+"/hapmaps/"+hapmap+"/"
logName    = hapmapDir+"process_log.txt"
inputFile  = hapmapDir+"haplotypeMap.txt"

with open(logName, "a") as myfile:
	myfile.write("*--------------------------------------------------*\n")
	myfile.write("| Log of 'hapmap.expand_definitions.py'            |\n")
	myfile.write("*--------------------------------------------------*\n\n")

#============================================================================================================
# Load haplotype map entries from 'haplotypeMap.txt' file.
#------------------------------------------------------------------------------------------------------------
# Initialize output file counter.
outCount = 0;
# Open 'haplotypeMap.txt'.
haplotypeMap_data = open(inputFile,'r')
# Process digested FASTA genome file, line by line.
while True:
	# haplotype map entries are sets of three lines.
	#    1) Parent dataset name.
	#    2) Child dataset name.
	#    3) LOH definitions.
	parent    = haplotypeMap_data.readline()	# C.albicans_SC5314
	parent    = parent[:-1]
	child     = haplotypeMap_data.readline()	# Ca_haploid
	child     = child[:-1]
	LOH_entry = haplotypeMap_data.readline()	# 2(b[1:495225] a[501338:3188548]); 3(a[1:1315157] b[1335699:2232035]); 4(a[1:844690] b[854160:1799406]); 5(a[1:1603444]);
	if not parent:								#     6(a[1:1190929]); 7(a[1:1033531]); 8(a[1:949617]); 9(a[1:2286390]);
		break  # EOF

	# Define output file for this haplotype entry.
	outputFile = hapmapDir+"haplotypeFragments."+str(outCount)+".txt"
	output     = open(outputFile, 'w')

	# trim last two characters off the LOH entry string.   ('; ')
	LOH_entry = LOH_entry[:-2]

	# break LOH string into per-chromosome definitions.
	LOH_chr_list = string.split(LOH_entry,"; ")
	for chr in LOH_chr_list:
		# trim last character off the chr entry string.   (')')
		chr       = chr[:-1]
		# split string into chrID and LOH_list.
		chr_parts = string.split(chr,"(")
		chrID     = chr_parts[0];
		LOH_list  = string.split(chr_parts[1]," ")
		for LOH_item in LOH_list:
			# trim last character off the string.   (']')
			LOH_item       = LOH_item[:-1]
			# split string into homologID and coordinates.
			LOH_item_parts = string.split(LOH_item,"[")
			homologID      = LOH_item_parts[0]
			coordinates    = string.split(LOH_item_parts[1],":")
			startbp        = coordinates[0]
			endbp          = coordinates[1]
			output.write(">"+child+".chr"+chrID+" ("+startbp+".."+endbp+") ("+str(int(endbp)-int(startbp)+1)+"bp) [*] "+homologID+"\n")
			output.write("NULL\n")

	# close output file.
	output.close();

	# increment output file counter.
	outCount += 1
