library(pbdMPI)
## Generate a set of three example protein sequences
seqs <- AAStringSet(c("seq1"="MLVVD",
                      "seq2"="PVVRA",
                      "seq3"="LVVR"))
## Count the kmers and generate a dataframe of the frequencies
freqs <- count_kmers(seqs, klen = 3, parallel = FALSE)
head(freqs)
##    kmer count
##1:  LVV  2
##2:  MLV  1
##3:  PVV  1
##4:  VRA  1
##5:  VVD  1
##6:  VVR  2
