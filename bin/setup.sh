#! /bin/bash

# Get the SRA toolkit if not in path
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.7.0/sratoolkit.2.7.0-centos_linux64.tar.gz
tar xzf sratoolkit.2.7.0-centos_linux64.tar.gz
export PATH=$PWD/sratoolkit.2.7.0-centos_linux64/bin:$PATH

# Download COG annotation

# Include repository directory in path



