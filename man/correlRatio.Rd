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

  V. Chabridon, L. Clouvel, B. Iooss, M. Il Idrissi and F. Robin, 2022,
  \emph{Variance-based importance measures in the linear regression context: Review, 
  new insights and applications}, Preprint
  
}

\author{
Bertrand Iooss
}

\examples{
x <- runif(100)
y <- round(x)
correlRatio(x,y)
}