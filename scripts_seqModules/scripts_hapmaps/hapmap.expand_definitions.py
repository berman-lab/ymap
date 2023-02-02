# Input Arguments
#   1) user     : 'darren'
#   2) hapmap   : 'test_Ca'
#   4) main_dir : '/home/bermanj/shared/links/'
#   5) logName  :
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

import string, sys;
user       = sys.argv[1];
hapmap     = sys.argv[2];
main_dir   = sys.argv[3];
hapmapDir  = main_dir+"users/"+user+"/hapmaps/"+hapmap+"/";
logName    = hapmapDir+"process_log.txt";
inputFile  = hapmapDir+"haplotypeMap.txt";

myfile = open(logName, "a");
myfile.write("*===============================================================================*\n");
myfile.write("| Log of 'scripts_seqModules/scripts_hapmaps/hapmap.expand_definitions.py'.     |\n");
myfile.write("*-------------------------------------------------------------------------------*\n");

#============================================================================================================
# Load haplotype map entries from 'haplotypeMap.txt' file.
#------------------------------------------------------------------------------------------------------------
# Initialize output file counter.
outCount = 0;
# Open 'haplotypeMap.txt'.
haplotypeMap_data = open(inputFile,'r');
# Process digested FASTA genome file, line by line.
while True:
	# haplotype map entries are sets of three lines.
	#    1) Parent dataset name.
	#    2) Child dataset name.
	#    3) LOH definitions.
	parent    = haplotypeMap_data.readline();
	parent    = parent.strip();
	child     = haplotypeMap_data.readline();
	child     = child.strip();
	LOH_entry = haplotypeMap_data.readline();
	if not parent:
		break;  # EOF

	myfile.write("| parent "+str(outCount)+":   "+parent+"\n");
	myfile.write("| child "+str(outCount)+":    "+child+"\n");
	myfile.write("| LOH entry "+str(outCount)+":"+LOH_entry+"\n");

	# Define output file for this haplotype entry.
	outputFile = hapmapDir+"haplotypeFragments."+str(outCount)+".txt";
	output     = open(outputFile, 'w');

	LOH_entry = LOH_entry.strip()	# trim off whitespace.
	LOH_entry = LOH_entry[:-1];	# trim last character (';') off the string.

	# break LOH string into per-chromosome definitions.
	LOH_chr_list = LOH_entry.split("; ");					# '3(b[1:2232036]); 9(a[1:1143195] b[1143196:2286390])' => ['3(b[1:2232036])', '9(a[1:1143195] b[1143196:2286390])']
	for chr in LOH_chr_list:
		chr       = chr[:-1];						# trim last character: '9(a[1:1143195] b[1143196:2286390])' => '9(a[1:1143195] b[1143196:2286390]'
		chr_parts = chr.split("(");					# split into chrID and LOH_list: '9(a[1:1143195] b[1143196:2286390]' => ['9', 'a[1:1143195] b[1143196:2286390]']
		chrID     = chr_parts[0];
		LOH_list  = chr_parts[1].split(" ");				# 'a[1:1143195] b[1143196:2286390]' => ['a[1:1143195]', 'b[1143196:2286390]']
		for LOH_item in LOH_list:
			LOH_item       = LOH_item[:-1];				# trim last character (']') off: 'b[1143196:2286390]' => 'b[1143196:2286390'
			LOH_item_parts = LOH_item.split("[");			# split into homologID and coordiantes: 'b[1143196:2286390' => ['b', '1143196:2286390']
			homologID      = LOH_item_parts[0];
			coordinates    = LOH_item_parts[1].split(":");		# split into coordinates: '1143196:2286390' => ['1143196', '2286390']
			startbp        = coordinates[0];
			endbp          = coordinates[1];
			output.write(">"+child+".chr"+chrID+" ("+str(startbp)+".."+str(endbp)+") ("+str(int(endbp)-int(startbp)+1)+"bp) [*] "+homologID+"\n");
			output.write("NULL\n");
			myfile.write("|\t>"+child+".chr"+chrID+" ("+str(startbp)+".."+str(endbp)+") ("+str(int(endbp)-int(startbp)+1)+"bp) [*] "+homologID+"\n");
			myfile.write("|\tNULL\n");

	output.close();	# close output file.
	outCount += 1	# increment output file counter.

myfile.write("*-------------------------------------------------------------------------------*\n");
myfile.write("| 'scripts_seqModules/scripts_hapmaps/hapmap.expand_definitions.py' completed.  |\n");
myfile.write("*===============================================================================*\n");
myfile.close();
