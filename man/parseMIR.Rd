% File man/parseMIR.Rd
\name{parseMIR}
\alias{parseMIR}
\title{Load the microRNA seeds from the mature microRNA file}
\description{
	Opens the mature microRNA files and extract the seed and the 1A seed of microRNAs
}
\usage{
	mirnas <- parseMIR("mature.miRNAs.txt")
}
\arguments{
	\item{file}{The name of the mature microRNAs file}
}
\seealso{
	\code{\link{anabel}}
}
\examples{
	file <- "mature.miRNAs.txt"
	mirnas <- parseMIR(file)
}
