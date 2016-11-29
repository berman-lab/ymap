# This script parallelizes mpileup calls by calling mpileup separately on each
# chromosome, and then concatenating the results.

import sys, os, time
import csv
import subprocess
from subprocess import Popen, call

samtools, fasta, bam, log_file, cores, genome_folder, final_pileup_file = sys.argv[1:]
cores = int(cores)

analysis_folder = os.path.split(bam)[0]

# In seconds, how much time to sleep before polling the samtools processes.
SLEEP = 1

# loading the relevant chromosomes from figure defenitions file
chroms = []
temp_pileup_files = []
with open(genome_folder + "figure_definitions.txt", "r") as figure_definitions_file:
    reader = csv.reader(figure_definitions_file, delimiter='\t')
    next(reader, None) # skip header line
    for line in reader:
    # the structure of line in the file:
    #	# Chr	Use	Label	Name	posX	posY	width	height
    # checking if Use equals 1  if yes adding chromosome name to list (equals original name in fasta
        if (line[1] == "1"):
            chroms.append(line[3])
            temp_pileup_files.append(os.path.join(analysis_folder,chroms[-1] + ".tmp.pileup"))

# Start the pileups going:
pending_chroms = list(chroms)
running_chroms = {}
while True:
    # Can we start running more pileups now?
    while len(running_chroms) < cores and len(pending_chroms) > 0:
        next_chrom = pending_chroms[0]
        del pending_chroms[0]
        pileup = \
            Popen([samtools, "mpileup",
                   "-f", fasta,
                   "-r", next_chrom, bam],
                  stdout=subprocess.PIPE,
                  close_fds=True)
        # As Ymap ignores quality scores, we throw out the last column.
        tmp_pileup_file = open(temp_pileup_files[chroms.index(next_chrom)], 'w')
        awk = \
            Popen(["awk", '{print $1 " " $2 " " $3 " " $4 " " $5}'],
                  stdin=pileup.stdout,
                  stdout=tmp_pileup_file,
                  close_fds=True)
        running_chroms[next_chrom] = (pileup, awk, tmp_pileup_file)

    # Wait for a little and check if any pileups finished:
    time.sleep(SLEEP)
    for chrom, (pileup, awk, tmp_file) in list(running_chroms.items()):
        if pileup.poll() is None or awk.poll() is None:
            continue
        pileup.stdout.close()
        tmp_file.close()
        del running_chroms[chrom]

    # Have we finished running all pileups?
    if len(running_chroms) == 0 and len(pending_chroms) == 0:
        break

# Concatenate the temp files into a single pileup:
call(["cat"] + temp_pileup_files,
     stdout=open(os.path.join(analysis_folder, final_pileup_file), 'w'))

# Remove the temporary files:
for temp_file in temp_pileup_files:
    os.remove(temp_file)
