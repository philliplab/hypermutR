Very very rough notes about implementation
==========================================

## Rough layout of the functions

.tallyMutInSeq
  in
    context
    enforce_context: ref query both
    context_must_match: T F
    invert_evolution: T F
    ref_seq
    query_seq
  out
    list(seq_name = list(
         num_mutated
         num_not_mutated
         mutated_pos
         not_mutated_pos)
  body
    loop along ref_seq
      check_contexts
      check_hypermut
      format_results

tallyHypermut
  in
    control_context
    hypermut_context
    query_seq -> query_alignment
    ref_seq -> NULLABLE
  out
    add to list:
      num_mut_hyper + vec of pos
      num_non_not_hyper + vec of pos
      same for control
      p-value
  body
    loop over contexts
      loop over seqs
        tally_mutations_in_seq
    loop over seq
      compute_fisher
    format_results

correctMut
  not supported until clarity about mutation direction

exclMut

dedupForHypermut

redupForHypermut

hypermutProcessFile
  read
  dedup
  tally
  ?correct : correct | exclude
  redup
  output

## Tests

The really important part is the tests.

A three way comparison is needed, comparing LANL vs Paul's implementation vs a
control which I understand perfectly.

A key fact to remember: It is a lot less effort to run a small number of big
files through LANL than a big number of small files - design the tests
accordingly.

Simulate the test sequences by chaining a bunch of patterns together seperated
by gaps.

Then run there controls through LANL and save the results:

list(name of imput data,
     md5sum of input data,
     data.frame(seq.num, 
                num.mut, 
                num.potential.mut, 
                num.control.mut,
                num.control.potential.mut))

Use seed values to ensure that the datasets are identical. Later store them as
.rda files.

Store the LANL results as dput code inside the test scripts and later migrate
them also to .rda files.


