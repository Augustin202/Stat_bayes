Assignment 3
================
Daniel Bonnery, Max Kramkimel, Augustin Poissonnier
Dec 5, 2023

# Bayesian Statistics, Assignment 3

### Install the following packages:

abind, dplyr, ggplot2, glmnet, invgamma, MASS, parallel, plyr, reshape2,
stats, stringr, tools, VGAM

### Copy and execute the project:

-   Clone the project to a local repository via:

``` r
system("git clone  --branch augustin https://github.com/Augustin202/Stat_bayes.git")
```

-   open the .rproj file in rstudio
-   run:

``` r
library(targets)
Sys.setenv(TAR_PROJECT = "assignment_3")
tar_make()
```

# Bayesian Statistics, Assignment 2

## Results

### Plot 1

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

We plot the histograms of the empirical means of q by values of s and
R2y. The red line represents s/k. The blue line represents the mean of
the empirical means. We use a log10 scale for the x axis.

### Plot 2

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

We plot the histograms of the empirical means of q by values of s and
R2y. The red line represents s/k. The blue line represents the mean of
the empirical means. We use a log10 scale for the x axis.

## How to run

### Install the following packages:

abind, dplyr, ggplot2, glmnet, invgamma, MASS, parallel, plyr, reshape2,
stats, stringr, SweaveLst, tools, VGAM

### Copy and execute the project:

-   Clone the project to a local repository via:

``` r
system("git clone https://github.com/Augustin202/Stat_bayes.git")
```

-   open the .rproj file in rstudio
-   run:

``` r
source(main.R)
```
