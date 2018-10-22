suppressMessages(library( "ape" )) # for "read.dna", "write.dna"
suppressMessages(library( "seqinr" )) # for "write.fasta" and "read.fasta"
# The conversion between ape and seqinr formats is tedious - I write and read to disk to do it.
suppressMessages(library( "hypermutR" ))

# TODO: REMOVE
# convenient debugging assignments
if (FALSE){
  fasta.file <- '/tmp/hd_seqs.fasta'
  seqinr::write.fasta(sequences = hd_seqs, names = names(hd_seqs), file.out = '/tmp/hd_seqs.fasta')
  output.dir <- '/tmp/hype'
  p.value.threshold <- 0.1
  fix.instead.of.remove <- TRUE
  fix.with <- 'R'
}

## Here is where the action is.
fasta.file <- Sys.getenv( "removeHypermutatedSequences_inputFilename" )
output.dir <- Sys.getenv( "removeHypermutatedSequences_outputDir" )
if( output.dir == "" ) {
    output.dir <- NULL
}
p.value.threshold <- Sys.getenv( "removeHypermutatedSequences_pValueThreshold" )
if( p.value.threshold == "" ) {
    p.value.threshold <- "0.1"
}
fix.instead.of.remove <- Sys.getenv( "removeHypermutatedSequences_fixInsteadOfRemove" );
if( ( nchar( fix.instead.of.remove ) == 0 ) || ( fix.instead.of.remove == "0" ) || ( toupper( fix.instead.of.remove ) == "F" ) || ( toupper( fix.instead.of.remove ) == "FALSE" ) ) {
    fix.instead.of.remove <- FALSE
} else {
    fix.instead.of.remove <- TRUE
}
fix.with <- Sys.getenv( "removeHypermutatedSequences_fixWith" )
if( fix.with == "" ) {
    fix.with <- "R"
}

# Generate names of output files with Paul's grep foo
if( length( grep( "^(.*?)\\/[^\\/]+$", fasta.file ) ) == 0 ) {
    fasta.file.path <- ".";
} else {
    fasta.file.path <-
        gsub( "^(.*?)\\/[^\\/]+$", "\\1", fasta.file );
}
fasta.file.short <-
    gsub( "^.*?\\/?([^\\/]+?)$", "\\1", fasta.file, perl = TRUE );
fasta.file.short.nosuffix <-
    gsub( "^([^\\.]+)(\\..+)?$", "\\1", fasta.file.short, perl = TRUE );
fasta.file.short.suffix <-
    gsub( "^([^\\.]+)(\\..+)?$", "\\2", fasta.file.short, perl = TRUE );
if( is.null( output.dir ) ) {
    output.dir = fasta.file.path;
}
if( is.null( output.dir ) ) {
    output.dir = ".";
}
## Remove "/" from end of output.dir
output.dir <-
    gsub( "^(.*?)\\/+$", "\\1", output.dir );

if( fix.instead.of.remove ) {
    out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_fixHypermutatedSequencesWith", toupper(fix.with), fasta.file.short.suffix, sep = "" );
} else {
    out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_removeHypermutatedSequences", fasta.file.short.suffix, sep = "" );
}

if( file.exists( fasta.file ) ) {
  in.fasta <- seqinr::read.fasta( fasta.file );

  x <- remove_hypermut(dat = in.fasta,
                       verbose = FALSE,
                       fix_with = ifelse(fix.instead.of.remove, fix.with, FALSE),
                       ancestor = 'consensus',
                       p_value = p.value.threshold)
  out.fasta <- x$seq_results

  # convert sequences to ape format
  tmp_file_name = tempfile()
  seqinr::write.fasta(sequences = out.fasta, names = names(out.fasta), file.out = tmp_file_name)
  result_in_ape_format <- ape::read.dna(tmp_file_name, format = 'fasta')

  ape::write.dna( result_in_ape_format, out.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
  write.csv(x$all_mut_pos, paste( output.dir, "/", fasta.file.short.nosuffix, "_potentialHypermutatedPositions.csv", sep = ''), row.names=F)

  print(length(x$seq_hypermutants))
} else {
  stop( paste( "File does not exist:", fasta.file ) );
}
