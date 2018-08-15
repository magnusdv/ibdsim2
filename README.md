
<!-- README.md is generated from README.Rmd. Please edit that file -->
ibdsim2 <img src="man/figures/logo.png" align="right" height=140/>
==================================================================

Introduction
------------

The purpose of ibdsim2 is to simulate and analyse the gene flow in pedigrees. In particular, such simulations can be used to study distributions of chromosomal segments shared *identical-by-descent* (IBD) by pedigree members. In each meiosis, the recombination process is simulated using sex specific recombination rates in the human genome (Kong et. al (2010) <doi:10.1038/nature09525>), or with recombination maps provided by the user. Simulations can be performed unconditionally, or by conditioning on specified IBD pattern. Additional functions provide summaries and further analysis of the allelic flow down through the pedigree.

ibdsim2 is an updated and improved version of [IBDsim](https://CRAN.R-project.org/package=IBDsim). In particular, the underlying pedigree structure is now imported from the [pedtools](https://github.com/magnusdv/pedtools) package instead of its predecessor [paramlink](https://CRAN.R-project.org/package=paramlink), which is no longer actively developed. In addition to the transition to pedtools, several new features are added in ibdsim2, including karyogram plots and analysis of *IBD absence* between (genealogically) related individuals.

Installation
------------

ibdsim2 is under development and not on CRAN yet. You can install the latest version from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install pedtools from github
devtools::install_github("magnusdv/ibdsim2")
```
