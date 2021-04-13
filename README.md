The following repository gathers the simulated datasets.

## Design 1

For the **Design 1** of simulations, data-sets can be easily simulated sourcing the script _simulate_desing_1.R_.
```r
  source("./simulate_desing_1.R")
```
It creates the function `get_design_1` with different transparent arguments and return a list with two named matrices, `X` ($n\times 1000$) and `Y` ($n\times 3$).


## Design 2

For the **Design 2** of simulations, data-sets can be easily simulated sourcing the script _simulate_desing_2.R_.
```r
  source("./simulate_desing_2.R")
```
In this script, the function `simulateMulti` is used as follows, where the argument `seed` fixes the seed ($\in [\![0,100]\!]$ for the scope of our simulations).
```r
seed <- 1
datas <- simulateMulti(seed=seed,plot=TRUE)
X <- cbind(Xs$X1,Xs$X2)
Y <- datas$Y
S_common <- cbind(datas$S$SXY$S1,datas$S$SXY$S2)
```
The option `plot`, set to `TRUE` in the example, allows to visualize the two blocks of the **x** part (first row) and the correlation between the **y** variable responses and the concatenation of the two blocks of the **x** part.
Plus, the `S_common` is the concatenation of the spectra, useful to get different statistics.
