#!/bin/sh

ml python/3.6.1
ml py-numpy
ml py-pandas
ml biology
ml py-biopython/1.79_py39

## Translate compactors and spindles to amino acid representation. 

seqkit translate -F --clean -f 6 $1 > $1"COMPACTORS.aa"

## Perform HMMer search of Pfam. 

Pfam_hmm=/oak/stanford/groups/horence/zheludev/preprint_060722/pipeline/Pfam/Pfam-A.hmm

hmmsearch --notextw -o PFAM.stdout \
    --tblout $1"_PFAM.tblout" --cpu 8 ${Pfam_hmm} $1"COMPACTORS.aa"
