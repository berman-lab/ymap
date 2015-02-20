#!/bin/bash -e
#
# project.ddRADseq.hapmap.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

user='darren1';
hapmap='C_albicans_SC5314_A21-s02-m09-r07';
project='test_ddRADseq_hapmap';

sh project.ddRADseq.hapmap.install_4.sh $user $project $hapmap;
