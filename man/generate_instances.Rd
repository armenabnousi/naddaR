% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nadda_functionized.R
\name{generate_instances}
\alias{generate_instances}
\title{Generate Instances for Indices in Protein Sequences}
\usage{
generate_instances(obj, labeled = TRUE, parallel = TRUE,
  nproc = ifelse(parallel, pbdMPI::comm.size(), 1), groundtruth = NULL,
  truth_filename = NULL, klen = 6, normalize = TRUE, impute = TRUE,
  winlen = 20, imputing_length = winlen\%/\%2, distributed = FALSE)
}
\arguments{
\item{obj}{A filepath to a fasta file containing protein sequences or an 
AAStringSet object containing the sequences}

\item{labeled}{TRUE if the method is called to construct a training set 
using a Pfam or InterPro labels or FALSE 
otherwise.}

\item{parallel}{Indicating whether the operation should be performed in 
parallel}

\item{nproc}{Currently not supported. Will use all processors available 
to the job on cluster}

\item{groundtruth}{A character string. Can be Pfam or InterPro}

\item{truth_filename}{The filepath to the labels file for generating the 
training set based on it}

\item{klen}{length of the k-mers to be used}

\item{normalize}{A boolean value, indicating whether the k-mer frequencies 
should be normalized}

\item{impute}{A boolean value, indicating whether imputed values should be 
inserted at the beginning and the 
end of the profiles}

\item{winlen}{An integer, size the window used for generation of each 
instance}

\item{imputing_length}{An integer, number of frequencies from the beginning 
and end of a sequence profile that should 
be used to impute the new values}

\item{distributed}{A boolean, indicating whether the data is spread among 
multiple processors.}
}
\value{
Returns a dataframe with one row for each instance. Each row 
contains \emph{winlen} k-mer frequencies around an 
index of a protein. The index number is stored in \emph{position} column. 
Name of the sequence is stored in \emph{name} 
column. 
If a training set is constructed, one column indicating whether it is a 
conserved index or not and a second column 
indicating the number of proteins in the dataset that have a similar 
conserved region are added to the returned 
dataframe.
}
\description{
Tconstructs a dataframe where each row corresponds to one 
index of one protein sequence from the 
input dataset. It can be used to generate training and test sets to train 
a NADDA classification model or to 
predict the conserved indices of input sequences based on a trained model.
}
\details{
Current version only supports Pfam and InterPro output files for 
generation of training set. The output 
from Pfam output file needs to be tabularized (replacing spaces with tabs).

If \emph{parallel} is set to \strong{TRUE} and  \emph{distributed} is set 
to \strong{FALSE}, the method distributes 
the data between different 
processors and sets \emph{distributed} to \strong{TRUE}. Otherwise, if the 
\emph{parallel} is set to \strong{FALSE} 
and \emph{distributed} is set to \strong{TRUE}, 
the kmer frequencies are computed on each processor separately but then 
communicated between each other, and therefore 
at the end all processors have the same set of frequencies for kmers 
stored, using which they will generate frequency 
profiles and instances of their chunk of sequences.
If you prefer to run the operation in serial, set both \emph{parallel} 
and \emph{distributed} to \strong{FALSE}.
}
\examples{
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

}
\author{
Armen Abnousi
}
