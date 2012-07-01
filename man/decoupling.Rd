\name{decoupling}
\alias{decoupling}
\alias{tell}
\alias{ask}

\title{Decoupling Simulations and Estimations}

\description{
  \code{tell} and \code{ask} are S3 generic methods for decoupling
  simulations and sensitivity measures estimations. In general, they are
  not used by the end-user for a simple \R model, but rather for an
  external computational code. Most of the sensitivity analyses objects
  of this package overload \code{tell}, whereas \code{ask} is overloaded
  for iterative methods only.
}

\usage{
tell(x, y = NULL, \dots)
ask(x, \dots)
}

\arguments{
  \item{x}{a typed list storing the state of the sensitivity study
  (parameters, data, estimates), as returned by sensitivity analyses
  objects constructors, such as \code{\link{src}}, \code{\link{morris}},
  etc.}
  \item{y}{a vector of model responses.}
  \item{\dots}{additional arguments, depending on the method used.}
}

\details{ 

  When a sensitivity analysis method is called with no model
  (i.e. argument \code{model = NULL}), it generates an incomplete object
  \code{x} that stores the design of experiments (field \code{X}),
  allowing the user to launch "by hand" the corresponding
  simulations. The method \code{tell} allows to pass these simulation
  results to the  incomplete object \code{x}, thereafter estimating the
  sensitivity measures.

%   For example, whith an external computational code, one
%   would typically:
%   \enumerate{
%     \item Save the data \code{x$X} to a file. Also, save the sensitivity
%     object \code{x} to a \code{.Rdata} file (see \code{\link{save}}).
%     \item Generate the simulator input files, for example with a Perl
%     script, or in \R with the  function \code{\link{template.replace}}.
%     \item launch the simulations, for example with a shell
%     script. Simulations may be distributed across a cluster of machines.
%     \item Gather the output files an extract the data of interest (shell
%     commands \code{grep} and \code{awk} are invaluable in
%     this task). We assume at this point that a file contains
%     a table of model responses.
%     \item Back to \R, load the sensitivity object \code{x} (saved at step 1) and
%     the responses \code{y} (saved at step 4).
%     \item Call the method \code{tell} to compute the sensitivity
%     measures, and do \code{print} and \code{plot} to analyze them.
%   }

  When the method is iterative, the data to simulate are not stored in
  the sensitivity analysis object \code{x}, but generated at each
  iteration with the \code{ask} method; see for example
  \code{\link{sb}}.

}

\value{
  \code{tell} doesn't return anything. It computes the sensitivity
  measures, and stores them in the list \code{x}.
  \strong{Side effect: \code{tell} modifies its argument \code{x}.}

  \code{ask} returns the set of data to simulate.
}

\keyword{methods}
