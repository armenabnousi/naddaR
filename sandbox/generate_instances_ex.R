library(Biostrings)
library(data.table)
## Generate a set of three example protein sequences
seqs <- AAStringSet(c("seq1"="MLVVD",
                      "seq2"="PVVRA",
                      "seq3"="LVVR"))
## Count the kmers and generate a dataframe of the frequencies
ins <- generate_instances(seqs, labeled = FALSE, parallel = FALSE, klen = 3, impute = TRUE, winlen = 5, normalize = FALSE)
head(ins)
##    name      position        freqn2  freqn1  freqindex       freqp1  freqp2
##    seq1      1               1.5     1.5     1.0             2.0     1.0
##    seq1      2               1.5     1.0     2.0             1.0     1.5
##    seq1      3               1.0     2.0     1.0             1.5     1.5
##    seq1      4               2.0     1.0     1.5             1.5     1.5
##    seq1      5               1.0     1.5     1.5             1.5     1.5
##    seq2      1               1.5     1.5     1.0             2.0     1.0

