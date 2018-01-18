% File man/revComp.Rd
\name{revComp}
\alias{revComp}
\title{Returns the reverse complement DNA sequence.}
\description{
	This function belongs to a R package called miRclump.
	It returns the reverse complementary of a DNA sequence.
}
\usage{
	revComp(seq="ATCGCTAGCTAGTC")
}
\arguments{
	\item{seq}{String of characters of DNA nucleotides.}
}
\seealso{
	\code{\link{miRclump}}
	\code{\link{parseMIR}}
}
\examples{
	revComp("ATGCGTGC")
	revComp(seq="ATGCGTGC")
	# [1] "GCACGCAT"
}
