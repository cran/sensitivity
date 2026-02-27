\name{correlRatio}
\alias{correlRatio}

\title{Correlation Ratio}

\description{
  \code{correlRatio} computes the correlation ratio between
  a quantitative variable and a qualitative variable
}

\usage{
correlRatio(X, y)
}

\arguments{
  \item{X}{a vector containing the quantitative variable.}
  \item{y}{a vector containing the qualitative variable (e.g. a factor).}
}

\value{
The value of the correlation ratio
}

\references{

  L. Clouvel, B. Iooss, V. Chabridon, M. Il Idrissi and F. Robin, 2024,
  \emph{An overview of variance-based importance measures in the linear regression context: 
  comparative analyses and numerical tests}, Socio-Environmental Systems Modelling, vol. 7, 
  18681, 2025, doi:10.18174/sesmo.1868.
  \url{https://hal.science/hal-04102053}
  
}

\author{
Bertrand Iooss
}

\examples{
x <- runif(100)
y <- round(x)
correlRatio(x,y)
}