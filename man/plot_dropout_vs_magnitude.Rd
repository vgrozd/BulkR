% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dropout_vs_magnitude.R
\name{plot_dropout_vs_magnitude}
\alias{plot_dropout_vs_magnitude}
\title{Title}
\usage{
plot_dropout_vs_magnitude(
  counts = NULL,
  expected_counts = NULL,
  dropout_fraq = NULL,
  per_cell = FALSE,
  logistic_fit = TRUE,
  fit.res = 300,
  color_by_libsize = TRUE,
  total_counts = NULL,
  force_sample = TRUE,
  seed = 142,
  ...
)
}
\arguments{
\item{counts}{A count vector/matrix to plot}

\item{expected_counts}{Vector of expected counts}

\item{dropout_fraq}{Array/vector of dropout fractions}

\item{per_cell}{Plot each single cell or average?}

\item{logistic_fit}{Fit a logistic curve?}

\item{fit.res}{Number of bins for the logistic curve}

\item{color_by_libsize}{Logical, plot line colors by lib size?}

\item{total_counts}{A vector of library sizes for each cell}

\item{force_sample}{Logical/numeric: Sample few cells to plot?}

\item{seed}{Numeric: Random sampling seed}

\item{...}{Graphic parameters for the plot lines}
}
\value{
Plot
}
\description{
Title
}
\examples{
plot_dropout_vs_magnitude(expected_counts(Seurat@assays$RNA@counts),
 dropout_fraq(
             Seurat@assays$RNA@counts,
   expected_counts(Seurat@assays$RNA@counts)),
   per_cell = FALSE,
   log = TRUE
 )
}
