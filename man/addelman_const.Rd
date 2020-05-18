\name{addelman_const}
\alias{addelman_const}
\title{Addelman and Kempthorne construction}
\description{
\code{addelman_const} implements the Addelman and Kempthorne construction of orthogonal arrays of strength two. 
}
\usage{
addelman_const(dimension, levels, choice="U")
}

\arguments{
  \item{dimension}{The number of columns of the orthogonal array.}
  \item{levels}{The number of levels of the orthogonal array. Either a prime number or a prime power number.}
  \item{choice}{A character from the list ("U","V","W","X") specifying which orthogonal array to construct (see "Details").}
  
}
\details{
The method of Addelman and Kempthorne allows to construct up to four orthogonal arrays. \code{choice} specify which orthogonal array is to be constructed. Note that the four orthognal arrays depends on each others through linear equations.
}

\value{
A matrix corresponding to the orthogonal array constructed.
}

\references{
A.S. Hedayat, N.J.A. Sloane and J. Stufken, 1999, \emph{Orthogonal Arrays: Theory and Applications}, Springer Series in Statistics.
}
\author{
Laurent Gilquin
}

\examples{
dimension <- 6
levels <- 7
OA <- addelman_const(dimension,levels,choice="U")
}
