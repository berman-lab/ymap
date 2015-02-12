#!/bin/bash -e
#
# project.ddRADseq.hapmap.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

user='default';
hapmap='C_albicans_SC5314_A21-s02-m09-r07';
project='Fig_06C.YJB12712-d2';

sh project.ddRADseq.hapmap.install_4.sh $user $hapmap $project;
