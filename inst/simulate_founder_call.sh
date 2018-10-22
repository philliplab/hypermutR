#!/bin/bash

echo $PWD

rm -f /tmp/hd_seqs.fasta || true
rm -rf /tmp/hype || true
mkdir /tmp/hype

R -q -e "suppressMessages(library(hypermutR)); seqinr::write.fasta(sequences = hd_seqs, names = names(hd_seqs), file.out = '/tmp/hd_seqs.fasta')"

export removeHypermutatedSequences_fixWith="R"
export removeHypermutatedSequences_fixInsteadOfRemove="0"
export removeHypermutatedSequences_pValueThreshold="0.1"
export removeHypermutatedSequences_inputFilename="/tmp/hd_seqs.fasta"
export removeHypermutatedSequences_outputDir="/tmp/hype"

R -q -f removeHypermutatedSequences.R --vanilla --slave

