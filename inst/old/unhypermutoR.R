library( "ade4", warn.conflicts = FALSE ) # needed by something.  ape?
library( "ape" ) # for "chronos", "as.DNAbin", "dist.dna", "read.dna", "write.dna"
library( "seqinr", warn.conflicts = FALSE ) # for "as.alignment", "consensus"

# This removes duplicate sequences from an input multiple alignment in fasta format, writes the output to another fasta file, writes a table of duplicate sequences, and returns the output filename.
# The output table is tab-delimited with two columns: retained and removed, where retained is a single sequence name (the one not removed) and removed is a comma-separated list of removed sequence names.
# increase.threshold.to.ensure.max = 100 will call maybe define as "duplicated" seqs that have HD <= k for some k > 0; the smallest k will be chosen such that the number of retained sequences is increase.threshold.to.ensure.max.
removeDuplicateSequencesFromAlignedFasta <- function ( input.fasta.file, output.dir = NULL, output.fasta.file = NULL, output.table.file = NULL, output.fasta.width = 72, add.copy.number.to.sequence.names = FALSE, include.gaps.in.Hamming = TRUE, increase.threshold.to.ensure.max = NULL ) {
# increase.threshold.until.at.most.n.seqs.remain
    if( length( grep( "^(.*?)\\/[^\\/]+$", input.fasta.file ) ) == 0 ) {
        input.fasta.file.path <- ".";
    } else {
        input.fasta.file.path <-
            gsub( "^(.*?)\\/[^\\/]+$", "\\1", input.fasta.file );
    }
    input.fasta.file.short <-
        gsub( "^.*?\\/?([^\\/]+?)$", "\\1", input.fasta.file, perl = TRUE );
    input.fasta.file.short.nosuffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\1", input.fasta.file.short, perl = TRUE );
    input.fasta.file.suffix <-
        gsub( "^([^\\.]+)(\\..+)?$", "\\2", input.fasta.file.short, perl = TRUE );

    if( !is.null( output.fasta.file ) ) {
      if( length( grep( "^(.*?)\\/[^\\/]+$", output.fasta.file ) ) == 0 ) {
          output.fasta.file.path <- NULL;
          output.fasta.file.path.is.absolute <- NA;
      } else {
          output.fasta.file.path <-
              gsub( "^(.*?)\\/[^\\/]+$", "\\1", output.fasta.file );
          output.fasta.file.path.is.absolute <- ( substring( output.fasta.file.path, 1, 1 ) == "/" );
      }
      output.fasta.file.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", output.fasta.file, perl = TRUE );

      if( !is.null( output.fasta.file.path ) && output.fasta.file.path.is.absolute ) {
          output.dir <- output.fasta.file.path;
      } else if( is.null( output.dir ) ) {
        if( is.null( output.fasta.file.path ) ) {
            output.dir <- input.fasta.file.path;
        } else {
            output.dir <- output.fasta.file.path;
        }
      } else {
          output.dir <- paste( output.dir, output.fasta.file.path, sep = "/" );
      }
      output.fasta.file <- output.fasta.file.short;
    } else { # is.null( output.fasta.file )
      output.fasta.file <- paste( input.fasta.file.short.nosuffix, "_removeDuplicateSequences", input.fasta.file.suffix, sep = "" );
    }
    if( is.null( output.dir ) || ( nchar( output.dir ) == 0 ) ) {
      output.dir <- ".";
    }
    ## Remove "/" from end of output.dir
    output.dir <-
      gsub( "^(.*?)\\/+$", "\\1", output.dir );

    if( is.null( output.table.file ) ) {
        output.table.file <- paste( input.fasta.file.short.nosuffix, "_removeDuplicateSequences.tbl", sep = "" );
    }
      if( length( grep( "^(.*?)\\/[^\\/]+$", output.table.file ) ) == 0 ) {
          output.table.file.path <- output.dir;
          output.table.file.path.is.absolute <- NA;
      } else {
          output.table.file.path <-
              gsub( "^(.*?)\\/[^\\/]+$", "\\1", output.table.file );
          output.table.file.path.is.absolute <- ( substring( output.table.file.path, 1, 1 ) == "/" );
      }
      output.table.file.short <-
          gsub( "^.*?\\/?([^\\/]+?)$", "\\1", output.table.file, perl = TRUE );
      output.table.file <- paste( output.table.file.path, output.table.file.short, sep = "/" );
    
    input.fasta <- read.dna( input.fasta.file, format = "fasta" );
    n.seqs <- nrow( input.fasta );

    # identical sequences have distance <= k between them [default k=0].
    # The pairwise.deletion = TRUE argument is necessary so that columns with any gaps are not removed.
    # The second call adds in the count of sites with gap differences
    # (gap in one sequence but not the other), which are not included
    # in the first call's "raw" count. Adding them back in, we are
    # effectively treating those as differences.

    in.dist <- dist.dna( input.fasta, model = "N", pairwise.deletion = TRUE );
    in.dist[ is.nan( in.dist ) ] <- 0;
    if( include.gaps.in.Hamming ) {
        # By default we include indels because we don't want to treat two sequences as identical when there are gaps between them.
        in.dist <- in.dist + dist.dna( input.fasta, model = "indel", pairwise.deletion = TRUE ); 
    }
    if( any( is.null( in.dist ) ) || any( is.na( in.dist ) ) || any( !is.finite( in.dist ) ) ) {
      ## TODO: REMOVE
      warning( "UH OH got illegal distance value" );
      print( "UH OH got illegal distance value" );
      print( in.dist );
    }
    dist.mat <- as.matrix( in.dist );
    ## NOTE: Because of ambiguity codes, the transitivity property
    ## might _not_ hold, eg you could have dist(A,B) == 0, dist(B,C) == 0, dist(A, C) == 1.

    # Ok so duplicates are those pairs with dist.mat <= threshold.k.
    # Keep the first.
    seq.removed.since.represented.by.seq <- rep( 0, nrow( dist.mat ) );
    names( seq.removed.since.represented.by.seq ) <- colnames( dist.mat );
    output.table.as.list <- list();
    removeDuplicatesAtThreshold <- function ( threshold.k ) {
      for( seq.i in ( 1:nrow( dist.mat ) ) ) {
          if( seq.removed.since.represented.by.seq[ seq.i ] > 0 ) {
              next;
          }
          # If the distances were true metric distances, ie if the transitivity property did hold, then it would be the case that we should not find this one to be identical to one already
          # removed, since then this one would also have already been
          # removed.  However, it does not hold (see above).  Thus we might need to add this to a collection of already removed sequences if it is zero-distance from something already removed.
          if( any( seq.removed.since.represented.by.seq != 0 ) && any( dist.mat[ seq.i, ( seq.removed.since.represented.by.seq > 0 ) ] <= threshold.k ) ) {
          ## TODO: REMOVE?
          #cat( "Found something within", threshold.k, "HD units from seq", seq.i, fill = TRUE );
              
              # This is part of a group that was already removed. It got
              # missed either because formerly threshold.k was lower, or perhaps even with constant threshold.k [eg this could happen if threshold.k == 0 by this being not perfectly identical to every other member due to ambiguities].  We go ahead and remove it if it is distance threshold.k from any member.
              .representative.seq <- seq.removed.since.represented.by.seq[ ( seq.removed.since.represented.by.seq > 0 ) & ( dist.mat[ seq.i, ] <= threshold.k ) ];
              stopifnot( length( .representative.seq ) > 0 );
              if( length( .representative.seq ) > 1 ) {
                  #if( threshold.k == 0 ) {
                      #stopifnot( all( .representative.seq == .representative.seq[ 1 ] ) );
                  #}
                  .representative.seq <- .representative.seq[ 1 ];
              }
              ## TODO: REMOVE
              # print( paste( "adding", colnames( dist.mat )[ seq.i ], "to", rownames( dist.mat )[ .representative.seq ] ) );
              output.table.as.list[[ rownames( dist.mat )[ .representative.seq ] ]][ 2 ] <-
                  paste( output.table.as.list[[ rownames( dist.mat )[ .representative.seq ] ]][ 2 ], colnames( dist.mat )[ seq.i ], sep = "," );
              seq.removed.since.represented.by.seq[ seq.i ] <<- .representative.seq;
          ## TODO: REMOVE
          #cat( "SETTING seq.removed.since.represented.by.seq[", seq.i, "] to", seq.removed.since.represented.by.seq[ seq.i ], fill = TRUE ); 
              next;
          }
          # Are there any "identical" to this one?
          if( !any( dist.mat[ seq.i, -seq.i ] <= threshold.k ) ) {
              next;
          }
          all.the.zeros <-
              colnames( dist.mat )[ -seq.i ][ dist.mat[ seq.i, -seq.i ] <= threshold.k ];
          output.table.as.list[[ rownames( dist.mat )[ seq.i ] ]] <<-
              c( colnames( dist.mat )[ seq.i ], paste( all.the.zeros, collapse = "," ) );
          #print( paste( all.the.zeros, collapse = ", " ) );
          ## TODO: REMOVE
          #cat( paste( "SETTING seq.removed.since.represented.by.seq[", all.the.zeros, "] to", seq.i, collapse = "\\n" ), fill = TRUE ); 
          seq.removed.since.represented.by.seq[ all.the.zeros ] <<- seq.i;
      } # End foreach seq.i

      return( NULL );
    } # removeDuplicatesAtThreshold ( threshold.k )
    if( is.null( increase.threshold.to.ensure.max ) ) {
        .result.ignored <- removeDuplicatesAtThreshold( 0 );
    } else {
      increase.threshold.to.ensure.max <-
          as.numeric( increase.threshold.to.ensure.max );
      stopifnot( increase.threshold.to.ensure.max > 1 );
      n.seqs.after <- n.seqs;
      threshold.k <- 0;
      ## TODO: REMOVE?
      cat( "TRYING TO ENSURE A MAXIMUM OF", increase.threshold.to.ensure.max, "sequences", fill = TRUE );
      cat( "THERE ARE", n.seqs.after, "sequences", fill = TRUE );
      #cat( n.seqs.after, ">", increase.threshold.to.ensure.max, "=", ( n.seqs.after > increase.threshold.to.ensure.max ), fill = TRUE );
      while( n.seqs.after > increase.threshold.to.ensure.max ) {
          n.seqs.before <- n.seqs.after;
          old.seq.removed.since.represented.by.seq <-
              seq.removed.since.represented.by.seq;
          old.output.table.as.list <- output.table.as.list;
          ## TODO: REMOVE?
          # cat( "n.seqs.before:", n.seqs.before, fill = TRUE );
          ## TODO: REMOVE?
          # cat( "removing dups at threshold:", threshold.k, fill = TRUE );
          # print( "sum( seq.removed.since.represented.by.seq )" );
          # print( sum( seq.removed.since.represented.by.seq ) );
          # print( "table( seq.removed.since.represented.by.seq )" );
          # print( table( seq.removed.since.represented.by.seq ) );
          # print( "sum( seq.removed.since.represented.by.seq == 0 )" );
          # print( sum( seq.removed.since.represented.by.seq == 0 ) );
          .result.ignored <- removeDuplicatesAtThreshold( threshold.k );
          
          #print( "sum( seq.removed.since.represented.by.seq )" );
          #print( sum( seq.removed.since.represented.by.seq ) );
          #print( "table( seq.removed.since.represented.by.seq )" );
          #print( table( seq.removed.since.represented.by.seq ) );
          #print( "table( seq.removed.since.represented.by.seq )" );
          #print( table( seq.removed.since.represented.by.seq ) );
          #print( "sum( seq.removed.since.represented.by.seq == 0 )" );
          #print( sum( seq.removed.since.represented.by.seq == 0 ) );
          # All the zeros, aka the number of unique remainders.
          n.seqs.after <- sum( seq.removed.since.represented.by.seq == 0 );
          ## TODO: REMOVE?
          #cat( "n.seqs.after:", n.seqs.after, fill = TRUE );
          if( n.seqs.after == 0 ) {
              ## TODO: REMOVE
              cat( "There are no sequences remaining at threshold", threshold.k, ", so we revert to", ( threshold.k - 1 ), ", which has", n.seqs.before, "sequences, despite the fact that it exceeds the limit (", increase.threshold.to.ensure.max, ").", fill = TRUE );
              ## woops. Undo it.
              threshold.k <- threshold.k - 1;
              n.seqs <- n.seqs.before;
              seq.removed.since.represented.by.seq <-
                  old.seq.removed.since.represented.by.seq;
              output.table.as.list <- old.output.table.as.list;
              break;
          } else if( n.seqs.after > increase.threshold.to.ensure.max ) {
              ## TODO: REMOVE
              cat( "There are still ", n.seqs.after, " sequences, which exceeds the limit (", increase.threshold.to.ensure.max, "). Increasing threshold to ", threshold.k + 1, ".", sep = "", fill = TRUE );
              if( n.seqs.after > increase.threshold.to.ensure.max ) {
                  threshold.k <- threshold.k + 1;
              }
          }
      } # End while( n.seqs.after > increase.threshold.to.ensure.max )
      # Consolodate seq.removed.since.represented.by.seq.
      # All the zeros, aka the number of unique remainders.
      n.seqs.after <- sum( seq.removed.since.represented.by.seq == 0 );
      ## TODO: REMOVE?
      cat( "DONE. Using threshold ", threshold.k, ": THERE ARE ", n.seqs.after, " sequences.", sep = "", fill = TRUE );
    } # End if is.null( increase.threshold.to.ensure.max ) .. else ..
    
    if( !is.null( output.table.file ) && ( length( output.table.as.list ) > 0 ) ) {
        output.table <- do.call( rbind, output.table.as.list );
        colnames( output.table ) <- c( "retained", "removed" );
        rownames( output.table ) <- NULL;
        write.table( output.table, file = output.table.file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE );
    }
    
    output.fasta <-
        input.fasta[ ( seq.removed.since.represented.by.seq == 0 ), , drop = FALSE ];
    if( add.copy.number.to.sequence.names ) {
        rownames( output.fasta ) <- sapply( rownames( output.fasta ), function ( .retained.seq.name ) {
            if( .retained.seq.name %in% names( output.table.as.list ) ) {
                #print( c( output.table.as.list[[ .retained.seq.name ]], gsub(  "[^,]", "", ( output.table.as.list[[ .retained.seq.name ]] )[ 2 ] ), nchar( gsub(  "[^,]", "", ( output.table.as.list[[ .retained.seq.name ]] )[ 2 ] ) ) + 2 ) );
                paste( .retained.seq.name, nchar( gsub(  "[^,]", "", ( output.table.as.list[[ .retained.seq.name ]] )[ 2 ] ) ) + 2, sep = "_" );
            } else {
                paste( .retained.seq.name, 1, sep = "_" );
            }
        } );
    }
               
    output.fasta.path <-
        paste( output.dir, "/", output.fasta.file, sep = "" );
    write.dna( output.fasta, output.fasta.path, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = output.fasta.width );

    # Return the file name.
    return( output.fasta.path );
} # removeDuplicateSequencesFromAlignedFasta ( input.fasta.file, ... )


#' Removes or corrects hypermutated sequences
#'
#' @param fasta.file The input sequnce file name.
#' @param output.dir Target output directory
#' @param p.value.threshold The cutoff for the Fisher Exact test to decide
#' whether or not hypermutation is present
#' @param fix.instead.of.remove If TRUE, then hypermutated sequences will not
#' be removed. The hyper mutated bases will be replaced with the value in
#' 'fix.with'.
#' @param fix.with The letter to replace hypermutated bases with. Default r.
#' @export

removeHypermutatedSequences <- function ( fasta.file, output.dir = NULL, p.value.threshold = 0.1, 
                                          fix.instead.of.remove = FALSE, fix.with = "r",
                                          consensus.file = NULL ) {

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

  if( fix.instead.of.remove ){
      fix.with <- tolower(fix.with)
  }
  
  in.fasta <- read.dna( fasta.file, format = "fasta" );

  if ( is.null(consensus.file) ){
      .consensus <- as.DNAbin( matrix( seqinr::consensus( as.character( in.fasta ) ), nrow = 1 ) );
      rownames( .consensus ) <- "Consensus sequence";
  } else {
      .consensus <- read.dna( consensus.file, format = "fasta" );
  }

  ## REMOVE DUPLICATES FIRST.
  ## If there are any duplicate sequences, remove them and
  ## incorporate the number of identical sequences into their names (_nnn suffix).
  # The output has the file name of the consensus file.
  fasta.file.no.duplicates <-
      removeDuplicateSequencesFromAlignedFasta( fasta.file, output.dir );
  fasta.file.no.duplicates.short <-
      gsub( "^.*?\\/?([^\\/]+?)$", "\\1", fasta.file.no.duplicates, perl = TRUE );
  fasta.file.no.duplicates.short.nosuffix <-
      gsub( "^(.*?)\\.[^\\.]+$", "\\1", fasta.file.no.duplicates.short, perl = TRUE );
  in.fasta.no.duplicates <-
      read.dna( fasta.file.no.duplicates, format = "fasta" );
  duplicate.sequences.tbl.file <-
      paste( output.dir, "/", fasta.file.no.duplicates.short.nosuffix, ".tbl", sep = "" );
  if( file.exists( duplicate.sequences.tbl.file ) ) { # won't exist if nothing was removed.
    duplicate.sequences.tbl.in <- read.table( file = duplicate.sequences.tbl.file, sep = "\t", header = TRUE );
    duplicate.sequences.tbl <- apply( duplicate.sequences.tbl.in, 1:2, function( .seq.name ) {
          # Fix our "escape" of the BAR and SLASH.
          .seq.name <- gsub( "-x-BAR-x-", "|", .seq.name );
          .seq.name <- gsub( "-x-SLASH-x-", "/", .seq.name );
          .seq.name <- gsub( "-x-BACKSLASH-x-", "\\", .seq.name );
          .seq.name <- gsub( "-x-DOT-x-", ".", .seq.name );
          .seq.name <- gsub( "-x-DOTDOT-x-", "..", .seq.name );
          return( .seq.name );
      } );
  } else {
    duplicate.sequences.tbl <- NULL;
  }
    
  compute.hypermut2.p.value <-
    function( seq.i, fix.sequence = FALSE ) {
      ## ERE I AM.  Bizarrely I can't quite reproduce what's on the web site.  I have more "potential" sites than are computed there, and I do not know why (todo: ask LANL)
      num.mut <- 0;
      num.potential.mut <- 0;
      num.control <- 0;
      num.potential.control <- 0;
      for( window.start.i in 1:( ncol( in.fasta.no.duplicates ) - 2 ) ) {
          # if the window has any gaps in either sequence, skip it.
          if( any( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0:2 ] ) == "-" ) || 
              any( as.character( .consensus[ 1, window.start.i + 0:2 ] ) == "-" ) ) {
              next;
          }
          #if( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ) == "a" ) && ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) ) {
          ## Added that the ref has to be an a or a g.
          if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) %in% c( "a", "g" ) ) && # 
              ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ) == "a" ) && 
              ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && 
              ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) && 
              ( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] == .consensus[ 1, window.start.i + 1 ] ) && 
              ( in.fasta.no.duplicates[ seq.i, window.start.i + 2 ] == .consensus[ 1, window.start.i + 2 ] ) ) { # unsure whether we need to enforce no change in the "context" sites.
              num.potential.mut <- num.potential.mut + 1;
              #if( ( as.character( .consensus[ window.start.i + 0 ] ) == "g" ) && ( as.character( .consensus[ window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( .consensus[ window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) ) {
              if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) == "g" ) ) { # don't enforce context in reference sequence
                  #print( window.start.i );
                  #print( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0:2 ] ) );
                  #print( as.character( .consensus[ 1, window.start.i + 0:2 ] ) );
                  num.mut <- num.mut + 1;
                  if( fix.sequence ) {
                      in.fasta.no.duplicates[ seq.i, window.start.i ] <<- as.DNAbin( fix.with );
                  }
              }
          }
          #if( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ) == "a" ) && ( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] ) %in% c( "c", "t" ) ) || ( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 2 ] ) == "c" ) ) ) ) {
          ## Added that the ref has to be an a or a g.
          if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) %in% c( "a", "g" ) ) && 
              ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0 ] ) == "a" ) && 
              ( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] ) %in% c( "c", "t" ) ) || 
                ( ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] ) %in% c( "a", "g" ) ) && 
                  ( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 2 ] ) == "c" ) ) ) && 
              ( in.fasta.no.duplicates[ seq.i, window.start.i + 1 ] == .consensus[ 1, window.start.i + 1 ] ) && 
              ( in.fasta.no.duplicates[ seq.i, window.start.i + 2 ] == .consensus[ 1, window.start.i + 2 ] ) ) { # unsure whether we need to enforce no change in the "context" sites.
              num.potential.control <- num.potential.control + 1;
              #if( ( as.character( .consensus[ window.start.i + 0 ] ) == "g" ) && ( as.character( .consensus[ window.start.i + 1 ] ) %in% c( "a", "g" ) ) && ( as.character( .consensus[ window.start.i + 2 ] ) %in% c( "a", "g", "t" ) ) ) {
              if( ( as.character( .consensus[ 1, window.start.i + 0 ] ) == "g" ) ) { # don't enforce context in reference sequence
                  #print( window.start.i );
                  #print( as.character( in.fasta.no.duplicates[ seq.i, window.start.i + 0:2 ] ) );
                  #print( as.character( .consensus[ 1, window.start.i + 0:2 ] ) );
                  num.control <- num.control + 1;
              }
          }
      }
      p.value <- fisher.test( ( matrix( c( num.control, ( num.potential.control - num.control ), num.mut, ( num.potential.mut - num.mut ) ), nrow = 2, byrow = T ) ) )$p.value;
      ## TODO: REMOVE
      print( c( num.mut = num.mut, num.potential.mut = num.potential.mut ) );
      print( c( num.control = num.control, num.potential.control = num.potential.control ) );
      return( p.value );
    }; # compute.hypermut2.p.value (..)
    
  exclude.sequence <- rep( FALSE, nrow( in.fasta ) );
  names( exclude.sequence ) <- rownames( in.fasta );
  fixed.sequence <- rep( FALSE, nrow( in.fasta ) );
  names( fixed.sequence ) <- rownames( in.fasta );
  for( seq.i in 1:nrow( in.fasta.no.duplicates ) ) {
    print(paste('Sequence: ', row.names(in.fasta.no.duplicates)[seq.i], sep = ''))
    p.value <- compute.hypermut2.p.value( seq.i );
    if( p.value < p.value.threshold ) {
        if( fix.instead.of.remove ) {
            .message <- paste( "Fixing", rownames( in.fasta.no.duplicates )[ seq.i ], "because the pseudo-HYPERMUT2.0 p-value is", p.value, "." );
            # Run it again but this time fix it.
            .result.ignored <-
                compute.hypermut2.p.value( seq.i, fix.sequence = TRUE );
            fixed.sequence[ rownames( in.fasta.no.duplicates )[ seq.i ] ] <- TRUE;
        } else {
            .message <- paste( "Excluding", rownames( in.fasta.no.duplicates )[ seq.i ], "because the pseudo-HYPERMUT2.0 p-value is", p.value, "." );
            exclude.sequence[ rownames( in.fasta.no.duplicates )[ seq.i ] ] <- TRUE;
        }
        ## TODO: REMOVE
        cat( .message, fill = TRUE );
    }
  } # End foreach seq.i

  if( is.null( duplicate.sequences.tbl ) ) {
    out.fasta <- in.fasta;
  } else {
    if( fix.instead.of.remove ) {
      out.fasta <- in.fasta;
    }
    #print( names( exclude.sequence )[ exclude.sequence ] );
    #print( duplicate.sequences.tbl[ , "retained" ] );
    excluded.sequences.with.duplicates <- intersect( names( exclude.sequence )[ exclude.sequence ], duplicate.sequences.tbl[ , "retained" ] );
    newly.excluded.sequence <- rep( FALSE, nrow( in.fasta ) );
    names( newly.excluded.sequence ) <- rownames( in.fasta );
    while( length( excluded.sequences.with.duplicates ) > 0 ) {
        newly.excluded.sequence[ newly.excluded.sequence ] <- FALSE;
        rownames( duplicate.sequences.tbl ) <- duplicate.sequences.tbl[ , "retained" ];
        .duplicates.strings <- duplicate.sequences.tbl[ excluded.sequences.with.duplicates, "removed" ];
        .duplicates.list <- strsplit( .duplicates.strings, "," );
        names( .duplicates.list ) <- excluded.sequences.with.duplicates;
        .result.ignored <- 
        lapply( excluded.sequences.with.duplicates, function( .retained.sequence ) {
            if( length( .duplicates.list[[ .retained.sequence ]] ) > 0 ) {
                if( fix.instead.of.remove ) {
                  ## TODO: REMOVE
                  cat( paste( "Fixing '", .duplicates.list[[ .retained.sequence ]], "' because it is a duplicate of ", .retained.sequence, ".", sep = "", collapse = "\n" ), fill = TRUE );
                  # Copy the fixed sequence.
                  .result.ignored <- 
                  sapply( .duplicates.list[[ .retained.sequence ]], function ( .retained ) {
                      out.fasta[ .retained, ] <<- out.fasta[ .retained.sequence, ];
                      fixed.sequence[ .duplicates.list[[ .retained.sequence ]] ] <<- TRUE;
                      return( NULL );
                  } );
              } else {
                  ## TODO: REMOVE
                  cat( paste( "Excluding '", .duplicates.list[[ .retained.sequence ]], "' because it is a duplicate of ", .retained.sequence, ".", sep = "", collapse = "\n" ), fill = TRUE );
              }
            }
            exclude.sequence[ .duplicates.list[[ .retained.sequence ]] ] <<- TRUE;
            newly.excluded.sequence[ .duplicates.list[[ .retained.sequence ]] ] <<- TRUE;
            return( NULL );
        } );
        # Recurse in case some of the sequences just excluded were recursively removed.
        excluded.sequences.with.duplicates <- intersect( names( newly.excluded.sequence )[ newly.excluded.sequence ], duplicate.sequences.tbl[ , "retained" ] );
    } # End while( length( excluded.sequences.with.duplicates ) > 0 )
    if( !fix.instead.of.remove ) {
        out.fasta <- in.fasta[ !exclude.sequence, , drop = FALSE ];
    }
  } # End if is.null( duplicate.sequences.tbl ) .. else ..
    
  # Write the subalignment as a fasta file
  if( fix.instead.of.remove ) {
      out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_fixHypermutatedSequencesWith", fix.with, fasta.file.short.suffix, sep = "" );
  } else {
      out.fasta.file = paste( output.dir, "/", fasta.file.short.nosuffix, "_removeHypermutatedSequences", fasta.file.short.suffix, sep = "" );
  }

  write.dna( out.fasta, out.fasta.file, format = "fasta", colsep = "", indent = 0, blocksep = 0, colw = 72 ); # TODO: DEHACKIFY MAGIC NUMBER 72 (fasta newline column)
    
  if( fix.instead.of.remove ) {
      .rv <- sum( fixed.sequence );
  } else {
      .rv <- sum( exclude.sequence );
  }
    
  return( .rv );
} # removeHypermutatedSequences (..)

