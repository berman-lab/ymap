#!/bin/bash

userProgramsLocation="/home/dabbey/software/";
#userProgramsLocation="~/software/";

bowtie2Directory=$userProgramsLocation"bowtie2-2.1.0/";
picardDirectory=$userProgramsLocation"picard-tools-1/picard-tools-1.105/";
fastqcDirectory=$userProgramsLocation"FastQC/";
gatkDirectory=$userProgramsLocation"gatk/GenomeAnalysisTK-2.8-1-g932cd3a/";
java7Directory=$userProgramsLocation"Java7/jdk1.7.0_51/jre/bin/";
seqtkDirectory=$userProgramsLocation"Seqtk/";

# Can be used to run PyPy (or any other Python implementation) instead of
# CPython for sripts that support it:
python_exec="/Applications/pypy-2.4.0-osx64/bin/pypy"