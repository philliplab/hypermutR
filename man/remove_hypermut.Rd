% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypermutR.R
\name{remove_hypermut}
\alias{remove_hypermut}
\title{Detect and process hypermutated sequences}
\usage{
remove_hypermut(dat, verbose = TRUE, fix_with = FALSE,
  ancestor = "consensus", p_value = 0.05)
}
\arguments{
\item{dat}{The sequence data. The structure must match the format produced
by read.fasta from the seqinr package. This is a list in which each element
represents a single sequence. Each element is of class \code{SeqFastadna}
and consists of a vector of single letters of class character with optional
attributes names and Annot.}

\item{verbose}{If TRUE, print the name and p-value of removed/fixed
sequences.}

\item{fix_with}{Either FALSE or a single letter. If not FALSE, then replace
the hypermutated base with the letter indicated.}

\item{ancestor}{Either 'consensus' to indicate that the consensus sequences
must be computed, or 'first' to indicate that the first sequence in the
dataset should be considered to be the ancestral sequence, or the ancestral
sequence itself.}

\item{p_value}{The p-value used by the one-sided fischer test.}
}
\value{
A list with three elements:
\describe{
  \item{all_mut_pos}{A data.frame that contains all the positions in all the
  sequences that are either a hypermutation or control position.}
  \item{seq_results}{A list that stores all the sequences that did not
  contain any hypermutation and, in the case that the fix_with parameter was
  set, those sequences with hypermutation that was corrected.}
  \item{seq_hypermutants}{A list of all the sequences that contains hypermutation.}
}
}
\description{
This function is a wrapper that will detect and either remove or 'fix'
hypermutated sequences depending of the value of the 'fix_with' argument.
}
\details{
It calls \code{ancestor_processing} to obtain the ancestral sequence to
compare the query sequences to, then calls \code{deduplicate_seqs} to remove
duplicate sequences for performance reasons, next loops over each unique
sequence, comparing it to the ancestral sequence with \code{scan_seq} and
finally collates the results.
}
\examples{
result <- remove_hypermut(hd_seqs)
names(result)
str(result)
}
\seealso{
\code{\link{ancestor_processing}}, \code{\link{deduplicate_seqs}},
and \code{\link{scan_seq}}
}
