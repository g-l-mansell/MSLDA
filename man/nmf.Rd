% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF_funcs.R
\name{nmf}
\alias{nmf}
\title{Run NMF using a count matrix}
\usage{
nmf(counts, K, max_iter = 50, thresh = 1e-04, seed = NULL)
}
\arguments{
\item{counts}{matrix of word counts}

\item{K}{internal dimension of matrix factors}

\item{max_iter}{the maximum number of iterations to run}

\item{thresh}{threshold for L convergence, (L_i - L_{i-1})/L_i < thresh}

\item{seed}{for the random initalisation of factors W and H}
}
\value{
A list of all parameters
}
\description{
Run NMF using a count matrix
}
