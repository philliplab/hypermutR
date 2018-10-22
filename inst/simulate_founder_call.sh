#!/bin/bash

rm /tmp/hd_seqs.fasta || true
rm -rf /tmp/hype || true
mkdir /tmp/hype

R -e "library(hypermutR); seqinr::write.fasta(sequences = hd_seqs, names = names(hd_seqs), file.out = '/tmp/hd_seqs.fasta')"

