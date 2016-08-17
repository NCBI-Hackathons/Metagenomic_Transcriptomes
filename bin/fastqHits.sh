#! /bin/bash

[[ -n "$1" ]] && sraid=$1 || exit
[[ -n "$2" ]] && cogid=$2 || exit

hitfile=$sraid.$cogid.out

export PATH=/home/ubuntu/sratoolkit.2.7.0-centos_linux64/bin:$PATH

spots=$(cat $hitfile | cut -f2 | cut -d'.' -f2 | sort -n | uniq | paste -d, -s -)

vdb-dump -f fastq1 -R $spots $sraid | \
    paste - - - - - - - - - - - - | \
    tee >(cut -f 1-4 | tr "\t" "\n" > ${sraid}.${cogid}_1.fq) | cut -f 9-12 | tr "\t" "\n" > ${sraid}.${cogid}_2.fq
