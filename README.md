# bstrl
Bayesian STreaming Record Linkage

This package performs streaming record linkage. We assume that files containing records about entities arrive sequentially in time. Each file is duplicate-free, 
but entities may be represented in more than one file. We want to determine, probabilistically, which records refer to the same entities across each file. We also want
these estimated links to be updated upon the arrival of each sequential file.

## Branches

* 'master' contains the latest version of streaming record linkage, using either Prior-Proposal-Recursive Bayes or Sequential Markov Chain Monte Carlo
* All other branches are dead ends or have been merged into 'master'

## Usage

To install, run

```
devtools::install_github("ianmtaylor1/bstrl", build_vignettes=TRUE)
```

then run `vignette()` to find included documentation and how-to's for this package.
