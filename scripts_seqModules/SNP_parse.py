# This script runs VariantsToTable for each file in the dataset vcf_files folder
# chromosome, and then concatenating the results.

import sys, os, time
import subprocess
from subprocess import Popen

gatk, fasta,  project_folder, log_file, cores = sys.argv[1:]
cores = int(cores)

# create subfolder to save the vcf's
vcf_folder = os.path.join(project_folder,"vcf_files")

filesToProcess = [] # will save the list of files to process

# get all files in the folder
for file in os.listdir(vcf_folder):
    if file.endswith(".vcf"):
        filesToProcess.append(os.path.join(vcf_folder,  file))
        
# In seconds, how much time to sleep before polling the samtools processes.
SLEEP = 1

# Start parsing:
pending_files = list(filesToProcess)
running_files = {}
tableFiles = [] # will save the name of the the table files created
while True:
    # Can we start running more pileups now?
    i = 0
    while len(running_files) < cores and len(pending_files) > 0:
        next_file = pending_files[0]
        tableFiles.append(os.path.splitext(filesToProcess[i])[0]  + '.table')
        del pending_files[0]
        parse = \
            Popen(["java", "-jar", gatk + "GenomeAnalysisTK.jar",
                   "-R", fasta,
                   "-T", "VariantsToTable",
                   "-V", filesToProcess[i] ,
                   "-F", "CHROM", 
                    "-F","POS", 
                    "-F", "REF", 
                    "-F",  "ALT", 
                    "-F", "AF",
                   "-o", os.path.join(tableFiles[i])],
                  stdout=subprocess.PIPE,
                  close_fds=True)
        running_files[next_file] = (parse)
        i += 1

    # Wait for a little and check if any pileups finished:
    time.sleep(SLEEP)
    for file, (parse) in list(running_files.items()):
        if parse.poll() is None:
            continue
        parse.stdout.close()
        del running_files[file]
    if len(running_files) == 0 and len(pending_files) == 0:
        break

# if we finished running VariantsToTable move on to concatenation
fout=open(os.path.join(vcf_folder,  "parsedVCF.txt"), "w+")
numOfFiles = len(tableFiles)
# first file:
for line in open(tableFiles[0]):
    fout.write(line)
# now the rest:    
for i in range(2,numOfFiles):
    f = open(tableFiles[i])
    f.next() # skip the header
    for line in f:
         fout.write(line)
    f.close() # not really needed
fout.close()

# remove all table files since they are not needed anymore
for tableFile in tableFiles:
    os.remove(tableFile)
