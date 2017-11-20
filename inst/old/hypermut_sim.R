library(ape)
library(devtools)
library(Biostrings)

ld_seqs <- readDNAStringSet('/home/phillipl/projects/hiv-founder-id/tests/sim_data/own_sim/HVTN503-162400146-1011_boost_1/HVTN503-162400146-1011_boost_1_own_sim.fasta')

hd_seqs <- readDNAStringSet('/home/phillipl/projects/hiv-founder-id/tests/sim_data/own_sim/HVTN503-162450071-1056_boost_2/HVTN503-162450071-1056_boost_2_own_sim.fasta')

devtools::use_data(ld_seqs)
devtools::use_data(hd_seqs)



# low d 1 highly mutated

n1 <- 1
n2 <- 0.9
n3 <- 0
ld1hm <- sim_hyper(ld_seqs, n1, n2, n3)
# h39c0_ATATGTATATGAATT_63

c_dat <- ld1hm
mutator <- which(grepl('^h[0-9]+c[0-9]+_', names(c_dat)))
names(c_dat[mutator])
as.character(c_dat[mutator])
as.character(ld_seqs[mutator])
as.character(hd_seqs[mutator])

writeXStringSet(ld1hm, '/home/phillipl/projects/hiv-founder-id/tests/sim_data/hypermut/ld1hm.fasta', width = 72)

# high d 1 highly mutated

n1 <- 1
n2 <- 0.9
n3 <- 0
hd1hm <- sim_hyper(hd_seqs, n1, n2, n3)
# h38c0_AATCCTCGAGTTGTC_12

writeXStringSet(hd1hm, '/home/phillipl/projects/hiv-founder-id/tests/sim_data/hypermut/hd1hm.fasta', width = 72)

# high d 1 30% seq

n1 <- 1
n2 <- 0.3
n3 <- 0
hd1hm_30 <- sim_hyper(hd_seqs, n1, n2, n3)
# h16c0_GGTTGAGTCGGGTTA_13

writeXStringSet(hd1hm_30, '/home/phillipl/projects/hiv-founder-id/tests/sim_data/hypermut/hd1hm_30.fasta', width = 72)

# high d 10 90% seq

n1 <- 10
n2 <- 0.9
n3 <- 0
hd1hm_90_10 <- sim_hyper(hd_seqs, n1, n2, n3)
# "h38c0_CTTCGGGGATTGTTT_7", "h38c0_ATATATTCAAGGCCA_4", "h38c0_AGGCAAGGGTCTTTT_9", "h38c0_GATCGGCTGTTGGGT_8", "h38c0_CATTAGCATAGCGAC_9", "h38c0_GGGAGCCTTTAGTTC_6", "h38c0_TGTAAGGTTCCCAGT_7", "h38c0_TTGGCACCCGTACTT_10", "h38c0_TAAGGGAGTCCTTTC_11", "h38c0_TTGGGCCTCTTGGCC_14"

writeXStringSet(hd1hm_90_10, '/home/phillipl/projects/hiv-founder-id/tests/sim_data/hypermut/hd1hm_90_10.fasta', width = 72)

# high d 1 70% control 100% seq

n1 <- 1
n2 <- 0.7
n3 <- 'all'
hd1hm_h70_c100_1 <- sim_hyper(hd_seqs, n1, n2, n3)
# h28c58_ATACGAGTGCTAGGA_8

writeXStringSet(hd1hm_h70_c100_1, '/home/phillipl/projects/hiv-founder-id/tests/sim_data/hypermut/hd1hm_h70_c100_1.fasta', width = 72)

# high d 1 all control all seq

n1 <- 1
n2 <- 'all'
n3 <- 'all'
hd1hm_h100_c100_1 <- sim_hyper(hd_seqs, n1, n2, n3)
# h43c54_CTGCGTAGTAGAATT_6

writeXStringSet(hd1hm_h70_c100_1, '/home/phillipl/projects/hiv-founder-id/tests/sim_data/hypermut/hd1hm_h100_c100_1.fasta', width = 72)

print('boo')
c_dat <- hd1hm_h100_c100_1
mutator <- which(grepl('^h[0-9]+c[0-9]+_', names(c_dat)))
names(c_dat[mutator])
as.character(c_dat[mutator])
as.character(ld_seqs[mutator])
as.character(hd_seqs[mutator])


n1 <- 1 # number of sequences to mutate
n2 <- 'all' # number of hypermutation sites to mutate
n3 <- 0 # number of control sites to mutate

dat <- ld_seqs
dat <- hd_seqs

test1 <- sim_hyper(dat, n1, n2='all', n3)
which(grepl('^h[0-9]+c[0-9]+_', names(test1)))
as.character(test1[710])
grepl('G[AG][AGT]', as.character(test1[710]))
grepl('G[CT]', as.character(test1[710]))
grepl('G[AG]C', as.character(test1[710]))

convert_alignment_to_matrix <- function(dat){
  stopifnot(length(unique(width(dat))) == 1) # all sequences of equal length
  mdat <- matrix('-', nrow = length(dat), ncol = width(dat)[1])
  for (i in 1:length(dat)){
    mdat[i,] <- strsplit(tolower(as.character(dat[i])), '')[[1]]
  }
  mdat
}

sim_hyper <- function(dat, n1, n2, n3){
  mdat <- convert_alignment_to_matrix(dat)
  if (length(n1) > 1){
    seqs_to_mut <- n1
  } else {
    seqs_to_mut <- sample(1:length(dat), n1)
  }
  for (i in seqs_to_mut){
    hypermut_pos <- numeric(0)
    control_pos <- numeric(0)
    for (j in 1:(ncol(mdat)-2)){
      if (mdat[i,j]   ==   'g' &
          mdat[i,j+1] %in% c('a', 'g') & 
          mdat[i,j+2] %in% c('a', 'g', 't')){
        hypermut_pos <- c(hypermut_pos, j)
      }
      if (mdat[i,j]   ==   'g' &
          mdat[i,j+1] %in% c('c', 't')){
        control_pos <- c(control_pos, j)
      }
      if (mdat[i,j]   ==   'g' &
          mdat[i,j+1] %in% c('a', 'g') & 
          mdat[i,j+2] ==   'c'){
        control_pos <- c(control_pos, j)
      }
    }
    # now do the n2 and n3 mutations
    if (n2 == 'all') {n2 <- length(hypermut_pos)}
    if (n3 == 'all') {n3 <- length(control_pos)}
    if (n2 < 1) {n2 <- floor(length(hypermut_pos)*n2)}
    if (n3 < 1) {n3 <- floor(length(control_pos)*n3)}
    if (length(hypermut_pos) < n2) {stop("not enough spots to hypermutate")}
    if (length(control_pos) < n3) {stop("not enough control spots to hypermutate")}
#    hypermut_pos_dist <- c(hypermut_pos_dist, length(hypermut_pos))
#    control_pos_dist <- c(control_pos_dist, length(control_pos))

    names(dat)[i] <- paste('h', n2, 'c', n3, '_', 
                           names(dat)[i], sep = '')
    print(names(dat)[i])

    n2_spots <- sample(hypermut_pos, n2)
    n3_spots <- sample(control_pos, n3)
    
    for (c_n2_spot in n2_spots){
      mdat[i, c_n2_spot] <- 'a'
    }
    for (c_n3_spot in n3_spots){
      mdat[i, c_n3_spot] <- 'a'
    }
  }
  new_dat <- DNAStringSet(apply(mdat, 1, paste, sep = '', collapse = ''))
  names(new_dat) <- names(dat)
  new_dat
}

