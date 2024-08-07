# mrnaScan

The Python toolkits designed to identify exons and 3'UTR of each reads for bulk/single-cell full-length mRNA-seq

## Installation
~~~
python3 -m pip install --upgrade mrnaScan
~~~

## Usage
~~~
# Basic usage
mrnaScan [options] -i <input.bam> -r <reference.gtf>

# For more information
mrnaScan -h
~~~

## Features
* The exons and 3'utr of each reads
* The pesudo-transcripts defined via the combination of exons
* The count table of pesudo-transcripts
* Other feature would be supported in the future ...

## Output
* *_exon_info.h5: a hdf5 file of count table of pesudo-transcripts and their exon information
* *_utr_info.tsv: a tab-separated table of 3'utr length of each gene
