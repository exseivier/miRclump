% File man/anabel.Rd
\name{anabel}
\alias{anabel}
\title{Generates anagrams for every microRNA seed}
\description{
	Generates anagrams for every microRNA seed including the 1A seed. Also it respects the 1A rule.
	This function takes a character vector which contains the microRNAs seeds and for every seed returns the reverse complement of the e seed, and 2:7 ading 1 adenine to the seed (1A rule).
}
\usage{
	miRanagrams <- anabel(mirnas)
}
\arguments{
	\item{mirnas}{'Character' vector which contains the microRNA seeds}
}
\seealso{
	\code{\link{parseMIR}}
	\code{\link{revComp}}
}
\examples{
	mirnas <- parseMIR("mature.miRNAs.txt")
	miRanagram <- anabel(mirnas)
}

