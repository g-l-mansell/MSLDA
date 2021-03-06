% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDA_reshaped.R, R/LDA_reshaped_noalpha.R
\name{lda_reshaped}
\alias{lda_reshaped}
\alias{lda_reshaped_noalpha}
\title{Run LDA adapted to use a count matrix}
\usage{
lda_reshaped(
  N,
  K,
  max_iter = 50,
  thresh = 1e-04,
  seed = NULL,
  cores = NULL,
  alpha = NULL
)

lda_reshaped_noalpha(
  N,
  K,
  max_iter = 50,
  thresh = 1e-04,
  seed = NULL,
  cores = NULL,
  alpha = NULL
)
}
\arguments{
\item{N}{matrix of word counts}

\item{K}{the number of topics to look for}

\item{max_iter}{the maximum number of EM iterations to run}

\item{thresh}{threshold for L convergence, (L_i - L_{i-1})/L_i < thresh}

\item{seed}{set a seed for the random documents to initialise beta}

\item{cores}{number of cores to run the E-step in parallel,
if NULL all detected cores are used}

\item{alpha}{if you want to set the exchangeable Dirichlet parameter for theta,
if NULL a default value of 1/K is used}
}
\value{
A list of all parameters

A list of all parameters
}
\description{
Run LDA adapted to use a count matrix
}
\section{Functions}{
\itemize{
\item \code{lda_reshaped_noalpha}: Alpha is fixed
}}

