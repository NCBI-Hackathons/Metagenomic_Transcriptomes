#! /bin/bash

[[ -n "$1" ]] && sraid=$1 || exit
[[ -n "$2" ]] && cogpath=$2 || exit

export PATH=/home/ubuntu/sratoolkit.2.7.0-centos_linux64/bin:$PATH

# Fetch the SRA data
prefetch $sraid

# Launch BLAST in parallel
ls $cogpath/COG????.fa | parallel --progress -j $(nproc) "tblastn_vdb -db $sraid -query {} -outfmt 6 -out $sraid.{/.}.out"

