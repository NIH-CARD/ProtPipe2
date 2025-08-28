R Package Guide: Installation
================

This guide is for R users who want to install the ProtPipe package and
use its functions in their own scripts.

## Prerequisites

- You need to have **R** installed on your system.
- The `devtools` package is required to install packages from GitHub.

If you do not have `devtools` installed, open R or RStudio and run the
following command:

``` r
install.packages("devtools")
```

## Installation Command

Once `devtools` is ready, you can install the ProtPipe package directly
from this GitHub repository with the following command:

``` r
devtools::install_github("NIH-CARD/ProtPipe2")
```

*(Note: Please replace the command above with the actual path to your
repository.)*

After the installation is complete, you can load the package into your R
session to start using it.

``` r
library(ProtPipe)
```

**Next Step: [Quick Start Example](./R-Package-Guide-2-Quick-Start)**
