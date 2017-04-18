# This script runs VariantsToTable for each file in the dataset vcf_files folder
# chromosome, and then concatenating the results.

import sys, os, time,  re
import subprocess
from subprocess import Popen

gatk, fasta, project_folder, log_file, cores = sys.argv[1:]
cores = int(cores)

# create subfolder to save the vcf's
vcf_folder = os.path.join(project_folder, "vcf_files")

files_to_process = [] # will save the list of files to process

# get all files in the folder
for file in os.listdir(vcf_folder):
    if file.endswith(".vcf"):
        files_to_process.append(os.path.join(vcf_folder, file))
        
# In seconds, how much time to sleep before polling the samtools processes.
SLEEP = 1

# Start parsing:
pending_files = list(files_to_process)
running_files = {}
table_files = [] # will save the name of the the table files created
i = 0
while True:
    # Can we start running more pileups now?
    while len(running_files) < cores and len(pending_files) > 0:
        next_file = pending_files[0]
        table_files.append(os.path.splitext(files_to_process[i])[0]  + '.table')
        del pending_files[0]
        parse = \
            Popen(["java", "-jar", gatk + "GenomeAnalysisTK.jar",
                   "-R", fasta,
                   "-T", "VariantsToTable",
                   "-V", files_to_process[i],
                   "-F", "CHROM", 
                   "-F", "POS", 
                   "-F", "REF", 
                   "-F", "ALT", 
                   "-F", "AF",
                   "-o", os.path.join(table_files[i])],
                   stdout=subprocess.PIPE,
                   close_fds=True)
        running_files[next_file] = parse
        i += 1

    # Wait for a little and check if any pileups finished:
    time.sleep(SLEEP)
    for file, parse in list(running_files.items()):
        if parse.poll() is None:
            continue
        parse.stdout.close()
        del running_files[file]
    if len(running_files) == 0 and len(pending_files) == 0:
        break

# if we finished running VariantsToTable move on to concatenation
with open(os.path.join(vcf_folder, "parsedVCF.txt"), "w+") as fout:
    num_of_files = len(table_files)
    # first file:
    first =  open(table_files[0])
    # write header
    fout.write(first.readline())
    for line in first:
        line_parts = re.split(r'\t+', line)
        # include only rows that contain SNPs
        if (len(line_parts[2]) == 1 and len(line_parts[3]) == 1):
            fout.write(line)
    # now the rest:    
    for i in range(1, numOfFiles):
        with open(table_files[i]) as f:
            f.next() # skip the header
            for line in f:
                line_parts = re.split(r'\t+', line)
                # include only rows that contain SNPs
                if (len(line_parts[2]) == 1 and len(line_parts[3]) == 1):
                    fout.write(line)
                    
# remove all table files since they are not needed anymore
for table_file in table_files:
    os.remove(table_file)
