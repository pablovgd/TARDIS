<!-- badges: start -->
[![R-CMD-check](https://github.com/pablovgd/T.A.R.D.I.S./actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pablovgd/T.A.R.D.I.S./actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# T.A.R.D.I.S. <img src="https://github.com/pablovgd/T.A.R.D.I.S./blob/main/www/tardis.png" width="150" height="150" align = right />        

R package for *TArgeted Raw Data Integration In Spectrometry*

## Installation
Make sure `R` (>=4.4.0) and `Rtools` are installed on your computer:

https://cloud.r-project.org/index.html

https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html

In R, run:

```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

devtools::install_github("pablovgd/T.A.R.D.I.S.",build_vignettes = TRUE,
  repos=c('http://cran.us.r-project.org',"http://bioconductor.org/packages/3.19/bioc"),
  dependencies=TRUE, type="source")

```

## Usage

Read the vignette for a tutorial on how to use `T.A.R.D.I.S.`

To launch the GUI in R:

```
library(TARDIS)
runTardis()

```

## Copyright

`T.A.R.D.I.S.` is licensed under the [GPLv3](http://choosealicense.com/licenses/gpl-3.0/)

As a summary, the GPLv3 license requires attribution, inclusion of copyright and license information, disclosure of source code and changes. Derivative work must be available under the same terms.

© Pablo Vangeenderhuysen (2024)
