#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -j y
#$ -v ARCH=LINUX
/nfs/home5/schles/md/run_MD-noPROFcon.pl fasta=3cat.f servers=yes 



