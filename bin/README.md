# __bin/__

## `inlineFastq2Fasta`

Converts fastq file on stdin to fasta file on stdout

## `blastCOGS.sh`

Uses GNU Parallel to blast each COG against an SRA accession.

```
usage: blastCOGS.sh [SRA accession] [path/to/cogs]
```

The `path/to/cogs` directory should contain one multi-fasta file for each functional group or COG. The sequences must be provided as amino acids. `blastCOGS.sh` will launch an instance of `tblastn_vdb` for each COG, with the protein sequences as the query and the SRA accession as the subject.

Example for setup on new Amazon node:

```bash
scp ubuntu@ip-172-16-242-19:/home/ubuntu/allcogs.tgz .
scp ubuntu@ip-172-16-242-19:/home/ubuntu/sratoolkit.2.7.0-centos_linux64.tar.gz .
git clone git@github.com:NCBI-Hackathons/Metagenomic_Transcriptomes.git

tar xzf sratoolkit.2.7.0-centos_linux64.tar.gz
tar xzf allcogs.tgz 

export PATH=$HOME/Metagenomic_Transcriptomes/bin:$PATH
export PATH=$HOME/sratoolkit.2.7.0-centos_linux64/bin:$PATH

. blastCOGS.sh SRR442395 $HOME/allcogs
```

## `fastqHits.sh`

Extract fastq from SRA for reads with hits.

```
usage: fastqHits.sh [SRA accession] [COG id]
```

The file `[SRAid].[COGid].out` is parsed to find spot numbers for reads with hits. The spot numbers are passed to `vdb-dump` and the output is split into `[SRAid].[COGid]_1.fq` and `[SRAid].[COGid]_2.fq`. We assume that there are 3 reads per spot and that the 2nd read is the barcode, need to implement other setups.

Example to extract reads and assemble using Trinity:

```bash
ls -S *.COG????.out | while read f; do
    v1=$(cut -d'.' -f1 <<<"$f")
    v2=$(cut -d'.' -f2 <<<"$f")
    fastqHits.sh $v1 $v2
    Trinity --seqType fq --max_memory 50G --left ${v1}.${v2}_1.fq --right ${v1}.${v2}_2.fq  --output ${v1}.${v2}.trinity &> ${v1}.${v2}.log
    [[ -e ${v1}.${v2}.trinity/Trinity.fasta ]] && \
        cp ${v1}.${v2}.trinity/Trinity.fasta ${v1}.${v2}.assembled.fasta || \
        touch ${v1}.${v2}.failed.
done
```

