#!/bin/bash -e
#
# project.ddRADseq.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

user='darren1';
project='test_ddRADseq_SC5314';

sh project.ddRADseq.install_4.sh $user $project;
