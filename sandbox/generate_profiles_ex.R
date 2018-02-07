library(Biostrings)
library(data.table)
## Generate a set of three example protein sequences
seqs <- AAStringSet(c("seq1"="MLVVD",
                      "seq2"="PVVRA",
                      "seq3"="LVVR"))
## Count the kmers and generate a dataframe of the frequencies
profs <- generate_profiles(seqs, klen = 3, parallel = FALSE, winlen = 5, normalize = FALSE)
head(profs)
profs
##[[1]]
##[[1]]$freqs
##[1] 1.5 1.5 1.0 2.0 1.0 1.5 1.5 1.5 1.5
##[[1]]$seq
##[1] "seq1"
##
##[[2]]
##[[2]]$freqs
##[1] 1.5 1.5 1.0 2.0 1.0 1.5 1.5 1.5 1.5
##[[2]]$seq
##[[1]] "seq2"
##
##[[3]]
##[[3]]$freqs
##[1] 2 2 2 2 2 2 2 2 
##[[3]]$seq
##[1] "seq3"
