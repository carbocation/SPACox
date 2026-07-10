# SPACox
A semiparametric empirical SPA approach based on a Cox regression model fitting (SPACox), that is scalable for a genome-wide single-variant survival analysis in large cohorts

### How to install and load this package

SPACox requires a Rust toolchain with `cargo` available on `PATH`. Rust 1.63
or newer is supported. The package compiles its native empirical-CGF backend
during installation.

```{r}      
library(devtools)  # author version: 2.1.0
install_github("WenjianBi/SPACox")
library(SPACox)
?SPACox  # manual of SPACox package
```
Current version is 0.1.2. For older version and version update information, plesase refer to OldVersions/

The Rust CGF backend is used by default and automatically uses the threads
available to the process. The thread count can be limited explicitly, and the
stable R implementation remains available for numerical comparisons:

```{r}
obj.null = SPACox_Null_Model(
  survival::Surv(time, event) ~ Cov1 + Cov2,
  data = Phen.mtx,
  pIDs = Phen.mtx$ID,
  gIDs = rownames(Geno.mtx),
  cgf.backend = "rust",
  cgf.threads = 8
)

obj.null.reference = SPACox_Null_Model(
  survival::Surv(time, event) ~ Cov1 + Cov2,
  data = Phen.mtx,
  pIDs = Phen.mtx$ID,
  gIDs = rownames(Geno.mtx),
  cgf.backend = "R"
)
```

Please do not hesitate to contact me (wenjianb@umich.edu) if you meet any problem. Suggestions or comments are also welcome.

### Reference

Wenjian Bi, Lars G. Fritsche, Bhramar Mukherjee, Sehee Kim, Seunggeun Lee, A Fast and Accurate Method for Genome-Wide Time-to-Event Data Analysis and Its Application to UK Biobank. American Journal of Human Genetics (2020), https://doi.org/10.1016/j.ajhg.2020.06.003

## coxphf
This directory contains some discussions about R package coxphf, an R package to conduct Firth's correction for Cox regression.
