#!/bin/bash -e
#
# project.ddRADseq.hapmap.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

user='darren1';
project='test_ddRADseq_11461_vs_SC5314_and_novelHapmap';
hapmap='test_hapmap_2';

sh project.ddRADseq.hapmap.install_4.sh $user $project $hapmap;
