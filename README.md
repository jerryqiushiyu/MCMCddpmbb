# MCMCddpmbb

R package for a Dynamic Dirichlet Process Mixture Beta-Binomial Model and Post-Processing Methods.

The MCMCddpmbb package is updated to be compatible with R version 4.0.2.

Before installing this package, you will need to install the latest version of Rtools, Rcpp, RcppArmadillo, reshape2, and ggplot2. 

There are two ways to install this package.

First, you can run the code below in R:

devtools::install_github("jerryqiushiyu/MCMCddpmbb")

Second, you can also follow the two steps below in command line.

1. Compile the R papckage into a tar ball with the code:
   R CMD build MCMCddpmbb

2. Install the tar ball with the code:
   R CMD INSTALL MCMCddpmbb_1.0.tar.gz
