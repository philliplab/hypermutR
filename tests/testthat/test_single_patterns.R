library(hypermutR)

set.seed(2)
non_mut_hyper <- simTestSeq(1, 0, 0, 0, F, 20)
non_mut_hyper_cons <- consensusMatrix(non_mut_hyper)
non_mut_hyper_cons <- DNAStringSet(paste(row.names(non_mut_hyper_cons)[apply(non_mut_hyper_cons, 2, which.max)], sep='', collapse=''))
names(non_mut_hyper_cons)[1] <- 'consensus'
non_mut_hyper_with_ref <- c(non_mut_hyper_cons, non_mut_hyper)

writeXStringSet(non_mut_hyper_with_ref, '/tmp/non_mut_hyper_1.fasta')

non_mut_hyper_with_ref_nogap <- DNAStringSet(gsub('-', '', non_mut_hyper_with_ref))
writeXStringSet(non_mut_hyper_with_ref_nogap, '/tmp/non_mut_hyper_1_nogap.fasta')
