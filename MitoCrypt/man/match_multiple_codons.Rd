% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_multiple_codons.r
\name{match_multiple_codons}
\alias{match_multiple_codons}
\title{Match multiple codons}
\usage{
match_multiple_codons(codons, Seq)
}
\arguments{
\item{codons}{character vector of searched codons.}

\item{Seq}{biostrings DNAstring object serving as template for search.}
}
\value{
vector of codon start locations, named with codon.
}
\description{
This functions is wrapper of matchPattern from biostrings package. Allows for multiple 3-nt long sequences matching.
}
\examples{

}
\keyword{codon}
\keyword{match}
