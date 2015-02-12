#!/bin/bash -e
#
# project.ddRADseq.install_4.sh
#
set -e;
## All created files will have permission 760
umask 007;

user='default'
project='Fig_04.SC5314';

sh project.ddRADseq.install_4.sh $user $project;
