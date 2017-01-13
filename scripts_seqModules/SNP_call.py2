# This script parallelizes CNV creation calls by calling HaplotypeCaller separately on each
# segment, and saving the result to vcf_files folder in the dataset folder

import sys, os, time
import csv
import subprocess
from subprocess import Popen

gatk, fasta, bam, log_file, cores = sys.argv[1:]
cores = int(cores)

project_folder = os.path.split(bam)[0]
# create subfolder to save the vcf's
vcf_folder = os.path.join(project_folder,"vcf_files")
if not os.path.exists(vcf_folder):
    os.makedirs(vcf_folder)

# In seconds, how much time to sleep before polling the gatk processes.
SLEEP = 1


chromSegments = [] # will save the name of the chromosome
startBP = []
endBP = []
ploidy = []
vcf_files = []
# loading chr segments file and running SNP caller for each segment
with open(os.path.join(project_folder, "chr_segments.txt"), "r") as chr_segments_file:
    reader = csv.reader(chr_segments_file, delimiter='\t')
    next(reader, None) # skip header line
    for line in reader:
        # the structure of line in the file (assuming contains only used chromosomes):
        #	#chr	startBP	endBP	ploidy
        chromSegments.append(line[0])
        startBP.append(line[1])
        endBP.append(line[2])
        ploidy.append(line[3])

# Start the vcf's going:
pending_segments = list(chromSegments)
running_segments = {}
while True:
    # Can we start running more pileups now?
    i = 0
    while len(running_segments) < cores and len(pending_segments) > 0:
        next_segment = pending_segments[0]
        del pending_segments[0]
        snpCall = \
            Popen(["java", "-jar", gatk + "GenomeAnalysisTK.jar",
                   "-R", fasta,
                   "-T", "HaplotypeCaller",
                   "-I", bam,
                   "-L", chromSegments[i] + ":" + startBP[i] + "-" + endBP[i],
                   "-ploidy",  ploidy[i], 
                   "-o", os.path.join(vcf_folder,chromSegments[i] + ":" + startBP[i] + "-" + endBP[i] + ".vcf")],
                  stdout=subprocess.PIPE,
                  close_fds=True)
        running_segments[next_segment] = (snpCall)
        i += 1

    # Wait for a little and check if any pileups finished:
    time.sleep(SLEEP)
    for chrom, (snpCall) in list(running_segments.items()):
        if snpCall.poll() is None:
            continue
        snpCall.stdout.close()
        del running_segments[chrom]

    # Have we finished running all haplotype calling?
    if len(running_segments) == 0 and len(pending_segments) == 0:
        break
