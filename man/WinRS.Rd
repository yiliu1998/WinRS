% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WinRS.R
\name{WinRS}
\alias{WinRS}
\title{Estimation and inference for the win ratio for two hierarchical endpoints using the S score estimator (nonparametric maximum likelihood estimator)}
\usage{
WinRS(
  df_a,
  df_b,
  cov.adjust = FALSE,
  covar.name = NULL,
  t.horizon,
  out1.name,
  delta1.name,
  out2.name,
  boot = TRUE,
  n_boot = 500,
  seed = 4399
)
}
\arguments{
\item{df_a}{data frame for data of treatment group a (treated group)}

\item{df_b}{data frame for data of treatment group a (treated group), the structure and names of variables should match df_a}

\item{cov.adjust}{whether to adjust covariate for missingness of the second endpoint Y2 (the quality-of-life outcome), default is FALSE}

\item{covar.name}{a vector for names of columns of covariates in df_a and df_b; if cov.adjust=TRUE, one must specify a set of covariates}

\item{t.horizon}{time horizon, a positive finite constant specified as the administrative end time of the study}

\item{out1.name}{column name of the primary survival outcome in both df_a and df_b}

\item{delta1.name}{column name of the primary event indicator in both df_a and df_b (1 indicates event observed, while 0 indicates right censored)}

\item{out2.name}{column name of the second outcome (non-survival) in both df_a and df_b}

\item{boot}{whether to conduct bootstrap for variance estimation, default is TRUE}

\item{n_boot}{number of bootstrap replicates if boot=TRUE, default is 500}

\item{seed}{seed for generating random numbers used in set.seed() function when splitting data; default is 4399}
}
\description{
Estimation and inference for the win ratio for two hierarchical endpoints using the S score estimator (nonparametric maximum likelihood estimator)
}
