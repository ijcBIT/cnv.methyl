% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMPextr2.R
\name{DMPextr}
\alias{DMPextr}
\title{Extract DMPs, annotation and methylation difference for each contrast}
\usage{
DMPextr(
  fit2,
  ContrastsDM,
  Pvalue,
  beta_normalized,
  mDiff,
  ann,
  writeOut = TRUE
)
}
\value{
data.table
}
\description{
For each contrast extract a set of DMPs and add gene annotation and methylation values
}
\author{
Angelika.Merkel & Izar.de.Villasante
}
