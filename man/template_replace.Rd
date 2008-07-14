\name{template.replace}
\alias{template.replace}

\title{Replace Values in a Template Text}

\description{
  \code{template.replace} replaces keys within special markups with
  values in a so-called template file. Pieces of \R code can be put into
  the markups of the template file, and are evaluated during the
  replacement.
}

\usage{
template.replace(text, replacement, eval = FALSE,
                 key.pattern = NULL, code.pattern = NULL)
}

\arguments{
  \item{text}{vector of character strings, the template text.}
  \item{replacement}{the list values to replace in \code{text}.}
  \item{eval}{boolean, \code{TRUE} if the code within
  \code{code.pattern} has to be evaluated, \code{FALSE} otherwise.}
  \item{key.pattern}{custom pattern for key replacement (see below)}
  \item{code.pattern}{custom pattern for code replacement (see below)}
}

\details{
  In most cases, a computational code reads its inputs from a text
  file. A template file is like an input file, but where some missing
  values, identified with generic keys, will be replaced by specific
  values.

  By default, the keys are enclosed into markups of the form \code{$(KEY)}.

  Code to be interpreted with \R can be put in the template text. Pieces
  of code must be enclosed into markups of the form
  \code{@\{CODE\}}. This is useful for example for formating the key
  values (see example). For interpreting the code, set \code{eval = TRUE}.

  Users can define custom patterns. These patterns must be
  perl-compatible regular expressions (see \code{\link{regexpr}}.
  The default ones are:
  \preformatted{key.pattern = "\\\\$\\\\(KEY\\\\)"
code.pattern = "@\\\\{CODE\\\\}"}
  Note that special characters have to
  be escaped both (one for perl, one for \R).
}

\value{
}

\examples{
txt <- c("Hello $(name)!", "$(a) + $(b) = @{$(a)+$(b)}",
         "pi = @{format(pi,digits=5)}")
replacement <- list(name = "world", a = 1, b = 2)
# 1. without code evaluation:
txt.rpl1 <- template.replace(txt, replacement)
print(txt.rpl1)
# 2. with code evalutation:
txt.rpl2 <- template.replace(txt, replacement, eval = TRUE)
print(txt.rpl2)
}

\keyword{IO}