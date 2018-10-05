hypermutR
=========

Hypermutation is a phenomena that affects HIV-1 introducing large numbers of
mutations into some sequences. It manifests in the datasets as sequences in
which large numbers of Guanine was mutated to Adenine, specifically when that
Guanine was surrounded by a particular pattern. The hypermut 2.0 tool available
from https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html is a
frequently used tool to detect and remove hypermutated sequences. We wrote a
new implementation of the hypermut 2.0 algorithm in R, which is available in
the hypermutR package on CRAN.

The hypermut algorithm compares each sequence in an alignment to some ancestral
sequence (usually approximated by the consensus sequence of the alignment),
tallying the frequency of specific mutations. Hypermutation occurs when a G
which is followed by an A or G (denoted by R in the IUPAC convention) and then
by an A, G or T (denoted by a D in the IUPAC convention) mutates to an A. More
compactly, when GRD become ARD, the mutation is flagged as possibly due to
hypermutation. In order to distinguish between true hypermutation and the
generally expected level of mutation, a baseline must be established. The
baseline is established by tallying G to A mutations when the G is followed
immediately by either a C or T (denoted by a Y in the IUPAC convention) or when
the G is followed by an A or a G (denoted by R in the IUPAC convention) and
then a C. More compactly, when GY becomes AY or GRC becomes ARC, the mutations
are tallied as the baseline mutation rate against which the potential
hypermutations must be compared. 

A one-sided Fisher’s exact test is used to compare the proportion of GRD
positions that became ARD positions to the proportion of GY or GRC positions
that became either AY or ARC positions. When the p-value of the test is smaller
than some threshold, with the default set to 0.1 as in (Abrahams et al., 2009),
then the individual sequence is flagged as a hypermutant and either the
sequence is removed from the dataset, or the mutated bases (the A’s followed by
RD) are replaced by an R to indicate that we are uncertain whether the mutation
was a random mutation or if it was the result of hypermutation.

In order to be a position of interest (either a control or hypermutation
position), all that is required is a G in the ancestral sequence. To classify
the position into either a hypermutation or control position, only the query
sequence is considered. If the two positions following the position that
contains the G in the ancestral sequence matches RD, then it is a hypermutation
position, else it is a control position. The two downstream positions in the
ancestral sequence are not considered. This implies the assumption that the two
downstream positions in the ancestral sequence mutates before the position of
interest.

## Installation Instructions for Ubuntu

Make sure you have a recent version of R. Follow
the instructions in the following link to set up the correct repositiory for apt:
http://stackoverflow.com/questions/10476713/how-to-upgrade-r-in-ubuntu. 

Make sure that both r-base and r-base-dev is installed
```{sh}
sudo apt-get install r-base r-base-dev
```

Next, install devtools' depedancies with apt-get:
```{sh}
sudo apt-get install libssl-dev libxml2-dev libcurl4-gnutls-dev
```

Then, from within R, install devtools:
```{r}
install.packages('devtools', repo = 'http://cran.rstudio.com/')
```

Finally, install hypermutR:
From a local file:

```{r}
library(devtools)
install_local('/path/to/file/hypermutR_x.y.z.tar.gz')
```

Please note that you must use install_local from devtools - install.packages
will not work. Change /path/to/file to the path to the installation file on
your computer and x.y.z to match the installation file you have.

Or using the bit_bucket repo:
```{r}
library(devtools)
install_bitbucket('hivdiversity/hypermutR', 
  auth_user = 'username', password = 'password')
```

Lastly, hypermutR includes a script that can be run from the commandline. You
need to put this script somewhere convenient ('/usr/bin' for example)
```{r}
file.symlink(from = file.path(find.package('hypermutR'), 'hypermutR.R'),
             to = '/usr/bin')
```

## Usage

### Within R

```{r}
library(hypermutR)
help('remove_hypermutation')
```

This will display the help for the main function in hypermutR.

### From the command line

```{sh}
hypermutR -h
```

or (depending on your installation):

```{sh}
hypermutR.R -h
```

This will display help for all the options and an example call to hypermutR.

More details can be found in the vignette, which can be viewed online at [hypermutR vignette](https://github.com/philliplab/hypermutR/blob/master/inst/doc/hypermut_vignette.pdf)
