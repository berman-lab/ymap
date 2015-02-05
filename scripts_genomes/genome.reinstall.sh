#!/bi
#n/bash -e
#
# genome.install_6.sh
#
set -e;
## All created files will have permission 760
umask 007;

user="darren1";
genome="test_fasta";

sh genome.install_6.sh user genome;
