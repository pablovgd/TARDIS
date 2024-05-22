<!-- badges: start -->
[![R-CMD-check](https://github.com/pablovgd/T.A.R.D.I.S./actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pablovgd/T.A.R.D.I.S./actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# T.A.R.D.I.S. <img src="https://github.com/pablovgd/T.A.R.D.I.S./blob/main/www/tardis.png" width="150" height="150" align = right />        

R package for *TArgeted Raw Data Integration In Spectrometry*

## Installation
Make sure `R` (**version >= 4.4.0**) is installed on your computer:

https://cloud.r-project.org/index.html

For **Windows** users, `Rtools` should be installed as well:

https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html

For **Mac** users, please install Xcode and a GNU Fortran compiler, see:
https://mac.r-project.org/tools/index.html
For Xcode, you can run this line in your Mac OS terminal:
```
sudo xcode-select --install
```

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

Load the package:

```
library(TARDIS)
```


**Read the vignette** for a tutorial on how to use `T.A.R.D.I.S.`

```
RShowDoc("gui_tutorial",package = "TARDIS")
```

To launch the GUI in R:

```
runTardis()
```

## Copyright

`T.A.R.D.I.S.` is licensed under the [GPLv3](http://choosealicense.com/licenses/gpl-3.0/)

As a summary, the GPLv3 license requires attribution, inclusion of copyright and license information, disclosure of source code and changes. Derivative work must be available under the same terms.

© Pablo Vangeenderhuysen (2024)
