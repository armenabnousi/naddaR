#library(pbdMPI)
#library(data.table)
#library(Biostrings)
#' @importFrom pbdMPI comm.rank
#' @importFrom pbdMPI comm.size
#' @importFrom pbdMPI allgather
#' @importFrom pbdMPI finalize
#' @import Biostrings
#' @importFrom Biostrings readAAStringSet
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom data.table rbindlist
#' @importFrom data.table as.data.table
#' @importFrom data.table setorder
#' @importFrom stats aggregate
#' @importFrom utils read.csv
#' @importFrom Rdpack reprompt

parallel_read_AA <- function(fasta_filename, 
                             nproc = comm.size()) {
        #print("in parallel_read_AA")
        seqs <- readAAStringSet(fasta_filename)
        seqs <- distribute_seqs(obj = seqs, nproc = nproc)
}
distribute_seqs <- function(obj, nproc = comm.size()) {
        #if (class(obj) == "character") {
        #        seqs <- distribute_seqs.character(obj = obj, nproc = nproc)
        #} else if (class(obj) == "AAStringSet") {
        #        seqs <- distribute_seqs.AAStringSet(obj = obj, nproc = nproc)
        #}
        UseMethod("distribute_seqs", obj)
        #return(seqs)
}
distribute_seqs.character <- function(obj, 
                                      nproc = comm.size()) {
	#print("indistribute_seqs.char")
        seqs <- readAAStringSet(obj)
        seqs <- distribute_seqs(obj = seqs, nproc = nproc)
}
distribute_seqs.AAStringSet <- function(obj, 
                                        nproc = comm.size()) {
	print("indistribute_seq.AA")
        rank <- as.integer(comm.rank())
        l <- length(obj)
        share <- ceiling(l / nproc)
        if ( (nproc - 1) * share  > l) {
                if (rank == 0) {
                        used_proc <- ceiling(l / share)
                        print (paste("Number of processors", nproc, 
                                "is unnecessary. Please rerun using only 
                                as many as", used_proc, "processors."))
                }
                finalize()
                quit()
        }
        if (rank < nproc) {
                procstart <- ((rank) * share) + 1
                procend <- min(procstart + share - 1, l)
                obj <- obj[procstart:procend]
                rm(procstart, procend)
        } else {
                obj <- c()
                #class(seqs) <- "AAStringSet"
        }
        rm(rank, share, l, nproc)
        return(obj)
}

generate_seqs_kmers <- function(seqs, klen) {
        #print("ingenerate_seq_kmer")
        kmers <- data.frame(kmer = c())
        kmers <- rbindlist(lapply(seqs, function(seq, klen){
                last <- nchar(as.character(seq)) - klen + 1
                seq_kmers <- sapply(1:last, function(start, seq, klen) {
                        a <- base::substr(seq, start, start + klen - 1)
                        a
                }, seq, klen)
                seq_kmers <- unique(seq_kmers)
                kmers <- data.frame(kmer=seq_kmers)
        }, klen))
        return(kmers)
}
kmerize <- function(obj, klen = 6, parallel = TRUE, 
                    nproc = ifelse(parallel, comm.size(), 1), 
                    distributed = FALSE) {
        #if (class(obj) == "character") {
        #        kmers <- kmerize.character(obj = obj, klen = klen, parallel = parallel, 
        #                          nproc = nproc, distributed = distributed)
        #} else if (class(obj) == "AAStringSet") {
        #        kmers <- kmerize.AAStringSet(obj = obj, klen = klen, parallel = parallel, 
        #                            nproc = nproc, distributed = distributed)
        #}
        UseMethod("kmerize", obj)        
        #return(kmers)
}
kmerize.character <- function(obj, klen = 6, parallel = TRUE, 
                              nproc = ifelse(parallel, comm.size(), 1), 
                              distributed = FALSE) {
        if (parallel) {
                seqs <- parallel_read_AA(obj, nproc)
                distributed <- TRUE
        } else {
                seqs <- readAAStringSet(obj)
        }
        kmers <- kmerize(obj = obj, klen = klen, 
                                     parallel = parallel, 
                                     nproc = nproc, 
                                     distributed = distributed)
}
kmerize.AAStringSet <- function(obj, klen = 6, parallel = TRUE, 
                                nproc = ifelse(parallel, comm.size(), 1), 
                                distributed = FALSE) {
        #print("inkmerize.AA")
        if (parallel & !distributed) {
                obj <- distribute_seqs(obj = obj, nproc = nproc)
        }
        kmers <- generate_seqs_kmers(obj, klen)
}

gather_kmers <- function(local_names, local_counts) {
        #print("ingather_kmers")
        gathered_freqs <- list()
        gathered_freqs <- allgather(local_counts)
        gathered_kmers <- list()
        gathered_kmers <- allgather(local_names)
        gathered_freqs <- unlist(gathered_freqs)
        gathered_kmers <- unlist(gathered_kmers)
        sum_received_freqs <- function(gathered_freqs, gathered_kmers) {
                freqs_df <- data.frame(kmer=gathered_kmers, 
                                       count=gathered_freqs)
                freqs_df <- data.table(freqs_df)
                setkey(freqs_df, "kmer")
                freqs <- aggregate(freqs_df$count, 
                                   list(kmers = as.factor(freqs_df$kmer)), 
                                   FUN = sum)
                rm(freqs_df)
                colnames(freqs) <- c("kmer", "count")
                freqs <- data.table(freqs)
                setkey(freqs, "kmer")
                freqs
        }
        freqs <- sum_received_freqs(gathered_freqs, gathered_kmers)
        rm(gathered_kmers, gathered_freqs)
        rm(sum_received_freqs)
        return(freqs)
}


#' Counting k-mers in the dataset.
#' @description counts the number of times each k-mer appears in the 
#' dataset and returns a dataframe indicating 
#' these counts.
#' Each k-mer in each sequence is counted
#' at most once, i.e., if there are multiple occureneces of a k-mer in one 
#' sequence, only one of them is counted.
#' @param obj A filepath to a fasta file containing protein sequences or an 
#' AAStringSet object containing the sequences
#' @param parallel Indicating whether the operation should be p
#' erformed in parallel
#' @param nproc Currently not supported. Will use all processors 
#' available to the job on cluster
#' @param klen length of the k-mers to be used
#' @param distributed A boolean, indicating whether the data is 
#' spread among multiple processors.
#' @details If \emph{parallel} is set to \strong{TRUE} and  
#' \emph{distributed} is set to \strong{FALSE}, the method 
#' distributes the data between different 
#' processors and sets \emph{distributed} to \strong{TRUE}. 
#' Otherwise, if the \emph{parallel} is set to \strong{FALSE} 
#' and \emph{distributed} is set to \strong{TRUE}, 
#' the kmer frequencies are computed on each processor separately 
#' but then communicated between each other, and therefore 
#' at the end all processors have the same set of frequencies 
#' for kmers stored, using which they will generate frequency 
#' profiles for their chunk of sequences.
#' If you prefer to run the operation in serial, set both 
#' \emph{parallel} and \emph{distributed} to \strong{FALSE}.
#' @return Returns a dataframe with two columns. Each row includes 
#' one k-mer and an integer indicating the number of 
#' times that k-mer appears in the input dataset. Each k-mer in 
#' each sequence is counted
#' at most once, i.e., if there are multiple occureneces of a 
#' k-mer in one sequence, only one of them is counted.
#!' @seealso \code{\link{generate_instances}} for generation of 
#!' training and test instances for each index in each 
#!' sequence.\cr
#!' \code{\link{generate_profiles} for generation of sequence k-mer 
#!' frequency profiles for each sequence}\cr
#!' \code{\link[pbdMPI]{comm.size}} for writing a distributed data object to a
#!' single file
#' @author Armen Abnousi
#' @example sandbox/count_kmers_ex.R
#!' @references{\insertRef{abnousi2016fast}{naddaR}}
#' @export count_kmers
count_kmers <- function(obj, klen = 6, parallel = TRUE, 
                        nproc = ifelse(parallel, comm.size(), 1), 
                        distributed = FALSE) {
        #print("incount_kmers")
        #if(class(obj) == "character") {
        #        kmers <- count_kmers.character(obj = obj, klen = klen, 
        #                              parallel = parallel, 
        #                              nproc = nproc, 
        #                              distributed = distributed)
        #} else if(class(obj) == "AAStringSet") {
        #        kmers <- count_kmers.AAStringSet(obj = obj, klen = klen, 
        #                                parallel = parallel, 
        #                                nproc = nproc, 
        #                                distributed = distributed)
        #}
        UseMethod("count_kmers", obj)
        #return(kmers)
}
#nah' @S3method count_kmers character
count_kmers.character <- function(obj, klen = 6, parallel = TRUE, 
                                  nproc = ifelse(parallel, comm.size(), 1), 
                                  distributed = FALSE) {
        #print("incount_kmers.char")
        if (parallel) {
                seqs <- parallel_read_AA(obj, nproc)
                distributed <- TRUE
        } else {
                seqs <- AAStringSet(readAAStringSet(obj))
        }
        freqs <- count_kmers(obj = obj, klen = klen, 
                                         parallel = parallel, 
                                         nproc = nproc, 
                                         distributed = distributed)
}
#nah' @S3method count_kmers AAStringSet
count_kmers.AAStringSet <- function(obj, klen = 6, parallel = TRUE, 
                                    nproc = ifelse(parallel, comm.size(), 1), 
                                    distributed = FALSE) {
        #print("incount_kmers.AA")
        warning(paste("number of processors,", nproc, ", is not supported in 
                      count_kmers function. Using all available processors."))
        if (parallel & !distributed) {
                obj <- distribute_seqs(obj = obj, nproc = nproc)
                distributed <- TRUE
        }
        kmers <- kmerize(obj = obj, klen = klen, 
                         parallel = parallel, nproc = nproc, 
                         distributed = distributed)
        freqs <- table(kmers)
        local_names <- names(freqs)
        local_counts <- as.integer(freqs)
        if (distributed) {
                freqs <- gather_kmers(local_names, local_counts)
                colnames(freqs) <- c("kmer", "count")
        } else {
                freqs <- data.table(kmer = local_names, 
                                    count = local_counts)
                setkey(freqs, "kmer")
        }
        rm(local_names, local_counts, kmers)
        return(freqs)
}

generate_seqs_freq_profiles <- function(d, klen, freqs, normalize, impute, 
                                        imputed_length = 0, 
                                        imputing_length = 0) {
        #print("ingenerate_seqs_freq_profiles")
        profiles <- lapply(names(d), function(seqname, klen, freqs, d, 
                                              imputed_length, imputing_length) {
                seq <- toString(d[seqname])
                last <- nchar(seq) - klen + 1
                seq_kmers <- sapply(1:last, function(start, seq, klen) {
                        a <- base::substr(seq, start, start + klen - 1)
                        a
                }, seq, klen)
                profile <- as.vector(freqs[as.character(seq_kmers)]$count)
                if (normalize) {
                        profile <- profile/max(profile)
                }
                if (impute) {
                        na_imputer_begin <- mean(profile[1 : imputing_length])
                        na_imputer_end <- 
                                mean(profile[(length(profile) - 
                                                      imputing_length + 1) : 
                                                     length(profile)])
                        profile <- c(rep(na_imputer_begin, imputed_length), 
                                     profile, rep(na_imputer_end, 
                                                  imputed_length + klen - 1))
                        rm(na_imputer_begin, na_imputer_end)
                }
                seq_profile = list(freqs = profile, seq = seqname)
                rm(profile, seq, seq_kmers, last)
                seq_profile
        }, klen, freqs, d, imputed_length, imputing_length)
        return(profiles)
}

#' Generate Protein k-mer frequency profiles
#' @description constructs a dataframe where each row corresponds to one 
#' index of one protein sequence from the 
#' input dataset. It can be used to generate training and test sets to train 
#' a NADDA classification model or to 
#' predict the conserved indices of input sequences based on a trained model.
#' @param obj A filepath to a fasta file containing protein sequences or an 
#' AAStringSet object containing the sequences
#' @param parallel Indicating whether the operation should be 
#' performed in parallel
#' @param nproc Currently not supported. Will use all processors 
#' available to the job on cluster
#' @param klen length of the k-mers to be used
#' @param normalize A boolean value, indicating whether the k-mer 
#' frequencies should be normalized
#' @param impute A boolean value, indicating whether imputed values 
#' should be inserted at the beginning and the 
#' end of the profiles
#' @param winlen An integer, size the window used for generation 
#' of each instance
#' @param imputing_length An integer, number of frequencies from the 
#' beginning and end of a sequence profile that should 
#' be used to impute the new values
#' @param distributed A boolean, indicating whether the data is 
#' spread among multiple processors.
#' @details If \emph{parallel} is set to \strong{TRUE} and  
#' \emph{distributed} is set to \strong{FALSE}, the method 
#' distributes the data between different 
#' processors and sets \emph{distributed} to \strong{TRUE}. 
#' Otherwise, if the \emph{parallel} is set to \strong{FALSE} 
#' and \emph{distributed} is set to \strong{TRUE}, 
#' the kmer frequencies are computed on each processor separately 
#' but then communicated between each other, and therefore 
#' at the end all processors have the same set of frequencies for 
#' kmers stored, using which they will generate frequency 
#' profiles for their chunk of sequences.
#' If you prefer to run the operation in serial, set both \emph{parallel} 
#' and \emph{distributed} to \strong{FALSE}.
#' @return Returns a list with one vector for each protein sequence in 
#' the dataset. A vector for sequence \emph{s} 
#' contains |s| - klen + 1
#' indices if \emph{impute} is set to \strong{FALSE} (where |s| is the 
#' length of the sequence). Otherwise it will 
#' include one index for each position in the sequence but also 
#' \emph{winlen \%\\\% 2} indices at 
#' the beginning and end of each sequence.
#!' @seealso \code{\link{generate_instances}} for generation of 
#!' training and test instances for each index in each sequence\cr
#!' \code{\link[pbdMPI]{comm.size}} for writing a distributed data object to a
#!' single file
#' @author Armen Abnousi
#' @example sandbox/generate_profiles_ex.R
#!' @references{\insertRef{abnousi2016fast}{naddaR}}
#' @export
generate_profiles <- function(obj, klen = 6, parallel = TRUE, 
                              nproc = ifelse(parallel, comm.size(), 1), 
                              normalize = TRUE, 
                              impute = TRUE, 
                              winlen = 20, imputing_length = winlen %/% 2, 
                              distributed = FALSE) {
        #print("ingenerate_profiles")
	#if (class(obj) == "character") {
	#        profiles <- generate_profiles.character(obj = obj, klen = klen, 
	#                                    parallel = parallel, 
	#                                    nproc = nproc, 
	#                                    normalize = normalize, 
	#                                    impute = impute, 
	#                                    winlen = winlen, 
	#                                    imputing_length = imputing_length, 
	#                                    distributed = distributed)
	#} else if (class(obj) == "AAStringSet") {
	#        profiles <- generate_profiles.AAStringSet(obj = obj, klen = klen, 
	#                                      parallel = parallel, 
	#                                      nproc = nproc, 
	#                                      normalize = normalize, 
	#                                      impute = impute, 
	#                                     winlen = winlen, 
	#                                      imputing_length = imputing_length, 
	#                                      distributed = distributed)
	#}
        UseMethod("generate_profiles", obj)
        #return(profiles)
}
#nah' @S3method generate_profiles character
generate_profiles.character <- function(obj, klen = 6, parallel = TRUE, 
                                        nproc = ifelse(parallel, comm.size(), 1), 
                                        normalize = TRUE, 
                                        impute = TRUE, 
                                        winlen = 20, 
                                        imputing_length = winlen %/% 2, 
                                        distributed = FALSE) {
        print("ingenerate_profiles.char")
        if (parallel & !distributed) {
                seqs <- parallel_read_AA(obj, nproc)
                distributed <- TRUE
        } else {
                seqs <- readAAStringSet(obj)
                distributed <- FALSE
        }
        profiles <- generate_profiles(obj = seqs, klen = klen, 
                                                  parallel = parallel, 
                                                  nproc = nproc, 
                                                  normalize = normalize, 
                                                  impute = impute, 
                                                  winlen = winlen, 
                                                  imputing_length = imputing_length, 
                                                  distributed = distributed)
}
#nah' @S3method generate_profiles AAStringSet
generate_profiles.AAStringSet <- function(obj, klen = 6, parallel = TRUE, 
                                          nproc = ifelse(parallel, 
                                                         comm.size(), 1), 
                                          normalize = TRUE, 
                                          impute = TRUE, 
                                          winlen = 20, 
                                          imputing_length = winlen %/% 2,
                                          distributed = FALSE) {
        #print("ingenerate_profiles.AA")
        imputed_length <- winlen %/% 2
        if (parallel & !distributed) {
                obj <- distribute_seqs(obj = obj, nproc = nproc)
                distributed <- TRUE
        }
        freqs <- count_kmers(obj = obj, klen = klen, parallel = parallel, 
                             nproc = nproc, distributed = distributed)
        profiles <- generate_seqs_freq_profiles(obj, klen, freqs, 
                                                normalize, impute, 
                                                imputed_length, 
                                                imputing_length)
        profiles
}

generate_empty_vectors_for_domains <- function(seq_lengths) {
        #print("ingenerate_empty_vectors_for_domains")
        doms <- lapply(as.character(seq_lengths$seqs), 
                       function(name, lengths) {
                len <- lengths[name]$lens
                ret <- rep(0, len)
                rm(len)
                list(name = name, dom = ret)
        }, seq_lengths)
        return(doms)
}
generate_domain_vectors <- 
        function(seq_lengths, groundtruth, truth_filename) {
        #print("ingenerate_domain_vectors")
        seq_dom_vectors <- generate_empty_vectors_for_domains(seq_lengths)
        if (groundtruth == "Pfam") {
                labels <- read.csv(truth_filename, header = FALSE, sep = "\t")
                dom_count <- table(labels$V6)
        } else if (groundtruth == "InterPro") {
                labels <- read.csv(truth_filename, header = FALSE, sep = "\t")
                dom_count <- table(labels$V5)
        }
        seq_domains <- lapply(seq_dom_vectors, 
                              function(dom, labels, groundtruth, dom_count){
                labels <- labels[labels$V1 == dom$name,]
                if (dim(labels)[1] > 0) {
                        if (groundtruth == "Pfam") {
                                for (row in 1:nrow(labels)) {
                                        dom$dom[labels[row, "V4"]:
                                                        labels[row, "V5"]] <- 
                                                dom_count[[(labels[row, "V6"])]]
                                }
                        } else if (groundtruth == "InterPro") {
                                for (row in 1:nrow(labels)) {
                                        dom$dom[labels[row, "V7"]:
                                                        labels[row, "V8"]] <- 
                                                dom_count[[(labels[row, "V5"])]]
                                }
                        }
                } else {
                        print(paste(dom$name, "does not have dim"))
                }
                dom
        }, labels, groundtruth, dom_count)
        rm(labels, dom_count, seq_dom_vectors)
        names(seq_domains) <- sapply(seq_domains, function(domain){domain$name})
        return(seq_domains)
}
extract_seqs_lengths <- function(d) {
        #print("inextract_seqs_lengths")
        lengths <- sapply(d, function(seq) {
                seq <- toString(seq)
                len <- nchar(seq)
                len
        })
        lengths <- data.frame(seqs=names(lengths), lens=lengths)
        lengths <- setorder(as.data.table(lengths))
        setkey(lengths, "seqs")
        return(lengths)
}
generate_instances_internal <- function(profiles, imputed_length, seq_lengths, 
                                        labeled = FALSE, 
                                        domain_vectors = NULL) {
        #print(ingenerate_instances_internal")
        if (is.null(domain_vectors)) domain_vectors <- 1
        instances <- rbindlist(lapply(profiles, function(profile, 
                                                         domain_vectors, 
                                                         imputed_length, 
                                                         seq_lengths) {
                name <- profile$seq
                if (labeled) {
                        domain_vector <- domain_vectors[[name]]
                } else {
                        domain_vector <- NULL
                }
                length <- seq_lengths[name]$lens
                indices <- seq((imputed_length + 1), length + imputed_length)
                seq_instances <- rbindlist(lapply(indices, 
                                                  function(i, profile, 
                                                           domain_vector, 
                                                           imputed_length){
                        index <- i - imputed_length
                        if (!is.null(domain_vector)) {
                                is_domain <- 
                                        ifelse (domain_vector$dom[index] > 0, 
                                                1, 0)
                                instance <- 
                                data.frame(name = profile$seq, 
                                   position = index, 
                                   label = is_domain, 
                                   dom_count <- 
                                     domain_vector$dom[index], 
                                     as.list(profile$freqs[(i - imputed_length)
                                                           :(i + imputed_length)
                                                           ]))
                        } else {
                                instance <- 
                                        data.frame(name = profile$seq, 
                                                   position = index, 
                                               as.list(profile$freqs[(i - 
                                                      imputed_length):
                                                      (i + imputed_length)]))
                        }
                }, profile, domain_vector, imputed_length))
                rm(name, domain_vector, length)
                seq_instances
        }, domain_vectors, imputed_length, seq_lengths))
        if (labeled) {
                colnames(instances) <- c("name", "position", "label", 
                                         "dom_count", 
                                         paste("freqn", imputed_length:1, 
                                               sep = ""), "freqindex", 
                                         paste("freqp", 1:imputed_length, 
                                               sep = ""))
        } else {
                colnames(instances) <- c("name", "position", 
                                         paste("freqn", imputed_length:1, 
                                               sep = ""), 
                                         "freqindex", 
                                         paste("freqp", 1:imputed_length, 
                                               sep = ""))
        }
        return(instances)
}

#' Generate Instances for Indices in Protein Sequences
#' @description Tconstructs a dataframe where each row corresponds to one 
#' index of one protein sequence from the 
#' input dataset. It can be used to generate training and test sets to train 
#' a NADDA classification model or to 
#' predict the conserved indices of input sequences based on a trained model.
#' @param obj A filepath to a fasta file containing protein sequences or an 
#' AAStringSet object containing the sequences
#' @param labeled TRUE if the method is called to construct a training set 
#' using a Pfam or InterPro labels or FALSE 
#' otherwise.
#' @param parallel Indicating whether the operation should be performed in 
#' parallel
#' @param nproc Currently not supported. Will use all processors available 
#' to the job on cluster
#' @param groundtruth A character string. Can be Pfam or InterPro
#' @param truth_filename The filepath to the labels file for generating the 
#' training set based on it
#' @param klen length of the k-mers to be used
#' @param normalize A boolean value, indicating whether the k-mer frequencies 
#' should be normalized
#' @param impute A boolean value, indicating whether imputed values should be 
#' inserted at the beginning and the 
#' end of the profiles
#' @param winlen An integer, size the window used for generation of each 
#' instance
#' @param imputing_length An integer, number of frequencies from the beginning 
#' and end of a sequence profile that should 
#' be used to impute the new values
#' @param distributed A boolean, indicating whether the data is spread among 
#' multiple processors.
#' @details Current version only supports Pfam and InterPro output files for 
#' generation of training set. The output 
#' from Pfam output file needs to be tabularized (replacing spaces with tabs).
#' 
#' If \emph{parallel} is set to \strong{TRUE} and  \emph{distributed} is set 
#' to \strong{FALSE}, the method distributes 
#' the data between different 
#' processors and sets \emph{distributed} to \strong{TRUE}. Otherwise, if the 
#' \emph{parallel} is set to \strong{FALSE} 
#' and \emph{distributed} is set to \strong{TRUE}, 
#' the kmer frequencies are computed on each processor separately but then 
#' communicated between each other, and therefore 
#' at the end all processors have the same set of frequencies for kmers 
#' stored, using which they will generate frequency 
#' profiles and instances of their chunk of sequences.
#' If you prefer to run the operation in serial, set both \emph{parallel} 
#' and \emph{distributed} to \strong{FALSE}.
#' @return Returns a dataframe with one row for each instance. Each row 
#' contains \emph{winlen} k-mer frequencies around an 
#' index of a protein. The index number is stored in \emph{position} column. 
#' Name of the sequence is stored in \emph{name} 
#' column. 
#' If a training set is constructed, one column indicating whether it is a 
#' conserved index or not and a second column 
#' indicating the number of proteins in the dataset that have a similar 
#' conserved region are added to the returned 
#' dataframe.
#!' @seealso \code{\link{generate_profiles}} for generation of 
#!' frequency profiles for each sequence\cr
#!' \code{\link[pbdMPI]{comm.size}} for writing a distributed data object to a
#!' single file
#' @author Armen Abnousi
#' @example sandbox/generate_instances_ex.R
#!' @references{\insertRef{abnousi2016fast}{naddaR}}
#' @export
generate_instances <- function(obj, labeled = TRUE, parallel = TRUE, 
                               nproc = ifelse(parallel, comm.size(), 1), 
                               groundtruth = NULL, truth_filename = NULL, 
                               klen = 6, normalize = TRUE, 
                               impute = TRUE, winlen = 20, 
                               imputing_length = winlen %/% 2, 
                               distributed = FALSE) {
        #print("ingenerate_instances")
        #if (class(obj) == "character") {
        #        instances <- generate_instances.character(obj = obj, 
        #                                                  labeled = labeled, 
        #                                     parallel = parallel, 
        #                                     nproc = nproc, 
        #                                     groundtruth = groundtruth, 
        #                                     truth_filename = truth_filename, 
        #                                     klen = klen, normalize = normalize, 
        #                                     impute = impute, winlen = winlen, 
        #                                     imputing_length = imputing_length, 
        #                                     distributed = distributed)
        #} else if (class(obj) == "AAStringSet") {
        #        instances <- generate_instances.AAStringSet(obj = obj, 
        #                                       labeled = labeled, 
        #                                       parallel = parallel, 
        #                                       nproc = nproc, 
        #                                       groundtruth = groundtruth, 
        #                                       truth_filename = truth_filename, 
        #                                       klen = klen, 
        #                                       normalize = normalize, 
        #                                       impute = impute, winlen = winlen, 
        #                                       imputing_length = imputing_length, 
        #                                       distributed = distributed)
        #}
        UseMethod("generate_instances", obj)
        #return(instances)
}
#nah' @S3method generate_instances character
generate_instances.character <- function(obj, labeled = TRUE, parallel = TRUE, 
                                         nproc = ifelse(parallel, 
                                                        comm.size(), 1), 
                                         groundtruth = NULL, 
                                         truth_filename = NULL, 
                                         klen = 6, normalize = TRUE, 
                                         impute = TRUE, winlen = 20, 
                                         imputing_length = winlen %/% 2, 
                                         distributed = FALSE) {
        print("ingenerate_instances.char")
        if (parallel & !distributed) {
                seqs <- parallel_read_AA(obj, nproc)
                distributed <- TRUE
        } else {
                seqs <- readAAStringSet(obj)
        }
        generate_instances(obj = seqs, labeled = labeled, 
                           parallel = parallel, 
                           nproc = nproc, 
                           groundtruth = groundtruth, 
                           truth_filename = truth_filename, 
                           klen = klen, normalize = normalize, 
                           impute = impute, winlen = winlen, 
                           imputing_length = imputing_length, 
                           distributed = distributed)
}
#nah' @S3method generate_instances AAStringSet
generate_instances.AAStringSet <- function(obj, labeled = TRUE, parallel = TRUE, 
                                           nproc = ifelse(parallel, 
                                                          comm.size(), 1), 
                                           groundtruth = NULL, 
                                           truth_filename = NULL, 
                                           klen = 6, normalize = TRUE, 
                                           impute = TRUE, winlen = 20, 
                                           imputing_length = winlen %/% 2, 
                                           distributed = FALSE) {
        #print("ingenerate_instances.AA")
        imputed_length <- winlen %/% 2
        if (parallel & !distributed) {
                warning("In this version of the package, nproc is not supported. 
                        All available processors will be used.")
                obj <- distribute_seqs(obj = obj, nproc = nproc)
                distributed <- TRUE
        }
        profiles <- generate_profiles(obj = obj, klen = klen, 
                                      parallel = parallel, 
                                      nproc = nproc, 
                                      normalize = normalize, 
                                      impute = impute, 
                                      winlen = winlen, 
                                      imputing_length = imputing_length, 
                                      distributed = distributed)
        seq_lengths <- extract_seqs_lengths(obj)
        domain_vectors <- NULL
        if (labeled) {
                if ((groundtruth != "Pfam" & groundtruth != "InterPro") | 
                    is.null(truth_filename)) {
                        stop("for labeled instances you need to 
                                provide the groundtruth type (Pfam or 
                             InterPro) and the filename for the labels")
                }
                domain_vectors <- generate_domain_vectors(seq_lengths, 
                                                          groundtruth, 
                                                          truth_filename)
        }
        instances <- generate_instances_internal(profiles, imputed_length, 
                                                 seq_lengths, labeled, 
                                                 domain_vectors)
}