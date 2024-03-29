% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supp_functions.R
\name{ATACplot}
\alias{ATACplot}
\title{Plot Corrlation plot for potential bias factors}
\usage{
ATACplot(cells, target_fragments, df, fc, project, cormethod)
}
\arguments{
\item{cells}{path to a list of cell barcodes}

\item{target_fragments}{path to a target fragments.tsv.gz file}

\item{df}{dataframe}

\item{fc}{Fold change matrix}

\item{project}{Name of the project}

\item{cormethod}{Method to use in the correlation calculation}
}
\value{
Creates a correlation plot and returns on target cnv data
}
\description{
This function creates a correlation plot for potential bias factors using
ATAC data
}
