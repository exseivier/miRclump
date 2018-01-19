% File man/freq.Rd
\name{freq}
\alias{freq}
\title{Calculates the frequenciesof observing any microRNA seed}
\description{
	Calculates the frequencies of observing any microRNA seed contained in a character vector
}
\usage{
	miRfreq <- freq(patterns, seqs, width)
}
\arguments{
	\item{patterns}{MicroRNA seeds character vector}
	\item{seqs}{DNAStringSet object storing the sequences of interest}
	\item{width}{Numeric argument. It is the width of the window for the seed frequency analysis}
}
\seealso{
	\code{\link{parseMIR}}
	\code{\link{anabel}}
}
\examples{
}
