#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option("--input_file", 
              help = "Path to input file and file name"),
  make_option("--output_file", 
              help = "Path and name of output file"),
  make_option("--p_value",
              default = 0.05,
              help = "The p-value to apply to the one-sided fischer test"),
  make_option("--ancestor",
              default = "consensus",
              help = "Either 'consensus' to indicate that the consensus sequences must be computed, or 'first' to indicate that the first sequence in the dataset should be considered to be the ancestral sequence, or the ancestral sequence itself."),
  make_option("--fix_with",
              default = NULL,
              help = "If omitted, hypermutants will be removed. If a single letter is specified, then hypermutants will be 'corrected' by replacing the hypermutated base with the specified letter.")
)

if (FALSE){
  # some debugging code that does not get run
  opt <- list(input_file = "/home/phillipl/projects/hypermutR/data/ld_seqs.fasta", ancestor = "consensus", output_file = '/tmp/bla.fasta')
}

opt <- parse_args(OptionParser(option_list = option_list,
  description = "Cast all hypermutants out into the cold!",
  epilogue = "Example Call:
hypermutR.R --input_file=/path/to/file.fasta --output_file=/path/to/out_file.fasta --p_value=0.1 --ancestor=first --fix_with=r

hypermutR.R --input_file=/home/phillipl/projects/hypermutR/data/ld_seqs.fasta --output_file=/tmp/ld_seqs_test.fasta --ancestor=consensus

File extensions must be specified in lower case."))

suppressPackageStartupMessages(library("hypermutR"))

if(is.null(opt$fix_with)){opt$fix_with <- FALSE}
if(!grepl('.fasta$', opt$input_file)){stop('input_file must end in .fasta (lowercase)')}
if(!grepl('.fasta$', opt$output_file)){stop('ouput_file must end in .fasta (lowercase)')}

dat <- read.fasta(opt$input_file)
cdat <- remove_hypermut(dat, ancestor = opt$ancestor, fix_with=opt$fix_with, p_value = opt$p_value)
names(cdat)
write.fasta(cdat$seq_result, names(cdat$seq_result), opt$output_file)
if (!is.null(cdat$seq_hypermutants)){
  write.fasta(cdat$seq_hypermutants, names(cdat$seq_hypermutants), 
              gsub(".fasta", "_hypermutants.fasta", opt$output_file))
}
write.csv(cdat$all_mut_pos, gsub('.fasta', '_mut_pos.csv', opt$output_file), row.names=F)

