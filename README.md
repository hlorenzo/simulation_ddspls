# Simulation designs

The following repository gathers the simulated data-sets used to test the data-driven sparse partial least squares (**ddsPLS**) methodology, accessible via the GitHub repository
GitHub for example **hlorenzo/ddsPLS2**
```r
  # install.packages("devtools")
  devtools::install_github("hlorenzo/ddsPLS2", build_vignettes = TRUE)
```

The RMarkdown **simulation_ddspls.Rmd** and its **HTML** version **simulation_ddspls.html** allow to get more insights on the ways to simulate following one of the three simulation structures. 

In a nutshell, the files **functions.R** gathers 3 functions to generate datasets associated with each of the 3 dataset structures:

* `get_toy_example` for **Toy Example** simulation structure,
* `get_design_1` for **Design 1** simulation structure,
* `get_design_2`  for **Design 2** simulation structure.

The **Design 2** simulation structure has already been used for publication purpose in

Cloarec, O. (2014), *Can we beat over‒fitting?*, Journal of Chemometrics, 28, pages 610– 614, **doi: 10.1002/cem.2602**
