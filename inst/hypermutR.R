#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option("--file_name", 
              help = "Path to input file and file name"),
  make_option("--output_name", 
              help = "Path and name of output file")
)

opt <- parse_args(OptionParser(option_list = option_list,
  description = "Cast all hypermutants out into the cold!",
  epilogue = "Example Call:
hypermutR.R --file_name=/path/to/file.fasta --output_name=/path/to/out_file.fasta"))

suppressPackageStartupMessages(library("hypermutR"))

dat <- readDNAStringSet(opt$file_name)
cdat <- remove_hypermut(dat)
writeXStringSet(cdat, opt$output_name, width=80)

