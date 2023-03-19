\name{testHSIC}
\alias{testHSIC}
\alias{print.testHSIC}
\alias{plot.testHSIC}




\title{Tests of Independence based on the Hilbert-Schmidt Independence Criterion (HSIC)}






\description{ 
\code{testHSIC} allows to test independence among all input-output pairs \eqn{(Xi, Y)} after a preliminary sensitivity analysis based on \acronym{HSIC} indices. \code{testHSIC} takes an object of class \code{sensiHSIC} (produced by a prior call to the function \code{sensiHSIC} that estimates \acronym{HSIC} indices) and it returns the estimated p-values after testing independence among all input-output pairs. For each input-output pair, having access to the p-value helps the user decide whether the null hypothesis \eqn{H0}: "\eqn{Xi} and \eqn{Y} are independent" must be accepted or rejected. If the kernels selected in \code{sensiHSIC} are all \bold{characteristic}, \eqn{H0} can be rewritten "\eqn{HSIC(Xi, Y)=0}" and this paves the way to several test procedures.

Depending on the sample size and the chosen test statistic (either a \bold{U-statistic} or a \bold{V-statistic}), there are up to four different methods to test \eqn{H0}. The \bold{asymptotic test} is recommended when the sample size \eqn{n} is around a few hundreds (or more). When \eqn{n} is smaller, a \bold{permutation-based test} must be considered instead. As a general rule, permutation-based tests can always be applied but a much heavier computational load is to be expected. However, if \acronym{HSIC} indices were initially estimated with V-statistics, the \bold{Gamma test} is a parametric method that offers an enticing tradeoff.
}





\usage{
testHSIC(sensi, test.method = "Asymptotic", B = 3000,
         seq.options = list(criterion = "screening", alpha = 0.05,
                            Bstart = 200, Bfinal = 5000, Bbatch = 100, Bconv = 200,
                            graph = TRUE) )
                            
\method{print}{testHSIC}(x, ...)

\method{plot}{testHSIC}(x, ylim = c(0, 1), err, ...)
}










\arguments{

  \item{sensi}{An object of class \code{"sensiHSIC"} which is produced by a prior call to the function \code{sensiHSIC}. In particular, \code{sensi} must contain objects named \code{"KX"} (3D-array filled with all input Gram matrices), \code{"KY"} (output Gram matrix), \code{"HSICXY"} (estimated \acronym{HSIC} indices) and \code{"estimator.type"} (either \code{"U-stat"} or \code{"V-stat"}). In addition, if \code{sensi} results from a conditional sensitivity analysis, \code{sensi} must also contain objects named \code{"cond"} (list of options describing the conditioning event) and \code{"weights"} (normalized conditioning weights).}
  
  \item{test.method}{A string specifying the numerical procedure used to estimate the p-values of the HSIC-based independence tests. Available procedure include \code{"Asymptotic"} (asymptotic test), \code{"Permutation"} (permutation-based test), \code{"Seq_Permutation"} (sequential permutation-based test) and \code{"Gamma"} (Gamma test). 
  
    \itemize{
      \item If \code{sensi} contains V-statistics, the asymptotic test (resp. the Gamma test) is recommended for large (resp. small) sample sizes. Otherwise, permutation-based tests can be used as well.
      \item If \code{sensi} contains U-statistics, the Gamma test must not be employed. The asymptotic test is recommended for large sample sizes. Otherwise, permutation-based tests can be used as well. 
    }
  
  }
  
  \item{B}{Number of random permutations carried out on the output samples before the non-parametric estimation of p-values. Only relevant if \code{test.method="Permutation"}.}
  
  \item{seq.options}{A list of options guiding the sequential procedure. Only relevant if \code{test.method="Seq_Permutaion"}.
  
    \itemize{
      \item \code{criterion} is a string specifying the stopping criterion. Available criteria include \code{"screening"} (permutations stop as soons as the estimated p-values have sufficiently converged so that they can be compared to the reference threshold \code{alpha}), \code{"ranking"} (permutations stop as soon as the estimated p-values have sufficiently converged so that they can be ranked) and \code{"both"} (permutations stop as soon as the two previous criteria are fulfilled).
      \item \code{alpha} is a scalar value (between \eqn{0} and \eqn{1}) specifying the type I error (probability of wrongly accepting \eqn{H0}). Only relevant if \code{criterion} is \code{"screening"} or \code{"both"}.
      \item \code{Bstart} is the initial number of random permutations before the first criterion check.
      \item \code{Bfinal} is the maximum number of random permutations.
      \item \code{Bbatch} is the number of permutations at each new iteration of the sequential procedure.
      \item \code{Bconv} is the number of permutations that is used to determine whether convergence has already occured or not. For \code{criterion="screening"}, convergence is assumed to be reached if the positions of the estimated p-values with respect to \code{alpha} no longer evolve after the \code{Bconv} latest permutations. For \code{criterion="ranking"}, convergence is assumed to be reached if the rankings of the estimated p-values no longer evolve after the \code{Bconv} latest permutations.
      \item \code{graph} is a boolean indicating whether the estimated p-values have to be plotted against the number of permutations.
    }
    
  }
  
  \item{x}{An object of class \code{"testHSIC"} storing the parameters and results of independence testing.}
  
  \item{ylim}{A vector of two values specifying the y-coordinate plotting limits.}
  
  \item{err}{A scalar value (between \eqn{0} and \eqn{1}) specifying the reference type I error. This value is used to plot a vertical line.}
  
  \item{...}{Additional options.}

}








\value{
  \code{testHSIC} returns a list of class \code{"testHSIC"}. It contains \code{test.method}, \code{B} (for the permutation-based test), \code{seq.options} (for the sequential permutation-based test) and the following objects:
  
  \item{call}{The matched call.}
  
  \item{pval}{The estimated p-values after testing independence for all input-output pairs.}
  
  \item{prop}{A vector of two strings.
    \itemize{
      \item The first string indicates if the chosen test method is asymptotic or non-asymptotic.
      \item The second string indicates if the chosen test method is parametric or non-parametric.
    }
  }
  
  \item{family}{Only if \code{test.method} is \code{"Asymptotic"} or \code{"Gamma"}. \cr
  A string indicating the parametric family used to estimate p-values.}
  
  \item{param}{Only if \code{test.method} is \code{"Asymptotic"} or \code{"Gamma"}. \cr 
  A \eqn{2}-column (resp. \eqn{3}-column) matrix containing the parameters of the Gamma (resp. Pearson III) distributions used to estimate p-values.}
  
  \item{Hperm}{Only if \code{test.method="Permutation"}. \cr
  A \eqn{B}-column matrix containing simulated values of the test statistics after randomly permuting the output samples. Each column in \code{Hperm} corresponds to one random permutation.}
  
  \item{paths}{Only if \code{test.method="Seq_Permutation"}. \cr
  A matrix containing all estimated p-values over the sequential test procedure. The \eqn{i}-th row provides all estimates of the \eqn{i}-th p-value as the number of permutations increases. If one row ends with a sequence of missing values \code{NA}, it means permutations were stopped earlier for this input variable. This can only happen if \code{test.method=screening}.}

}










\details{

\subsection{Why and how to properly choose kernels?}{

For a given input-output pair of variables, the Hilbert-Schmidt independence criterion \acronym{(HSIC)} is a dissimilarity measure between the joint bivariate distribution and the product of marginal distributions. Dissimilarity between those two distributions is measured through the squared norm of the distance between their respective embeddings in a reproducing kernel Hilbert space (\acronym{RKHS}) that directly depends on the selected input kernel \eqn{Ki} and the selected output kernel \eqn{KY}.

It must always be kept in mind that this criterion allows to detect independence within the pair \eqn{(Xi, Y)} provided that the two kernels are \bold{characteristic}.  
  \itemize{
    \item If both kernels are characteristic, \eqn{H0}: "\eqn{Xi} and \eqn{Y} are independent" is equivalent to \eqn{H0}: "\eqn{HSIC(Xi, Y)=0}" and any estimator of \eqn{HSIC(Xi, Y)} emerges as a relevant test statistic.
    \item If they are not, testing \eqn{H0}: "\eqn{HSIC(Xi, Y)=0}" is no longer sufficient for testing  \eqn{H0}: "\eqn{Xi} and \eqn{Y} are independent".
  }
The reader is referred to Fukumizu et al. (2004) for the mathematical definition of a characteristic kernel and to Sriperumbur et al. (2010) for an overview of the major related results. 

Responsability for kernel selection is left to the user while calling the function \code{sensiHSIC}. Let us simply recall that:
  \itemize{
    \item The Gaussian kernel, the exponential kernel, the Matern \eqn{3/2} kernel and the Matern \eqn{5/2} kernel (all defined on \eqn{R^2}) are \bold{characteristic}. They remain \bold{characteristic} when they are restricted to a compact domain \eqn{D} within \eqn{R^2}.
    \item The transformed versions of the four abovementioned kernels (all defined on \eqn{[0,1]^2}) are \bold{characteristic}.
    \item All Sobolev kernels (defined on \eqn{[0,1]^2}) are \bold{characteristic}.
    \item The categorical kernel (defined on any discrete probability space) is \bold{characteristic}.
  }

}

\subsection{Which test method is most appropriate?}{

The test statistic for the pair \eqn{(Xi, Y)} is either the \bold{U-statistic} or the \bold{V-statistic} associated to \eqn{HSIC(Xi, Y)}.

If a \bold{V-statistic} was used in \code{sensiHSIC}, four different test methods can be considered.
  
  \itemize{
  
    \item The \bold{asymptotic test} can be used if the sample size \eqn{n} is large enough (at least a hundred of samples). The asymptotic distribution of the test statistic is approximated by a Gamma distribution whose parameters are estimated with the method of moments. See Gretton et al. (2007) for more details about how to estimate the first two moments of the asymptotic Gamma distribution.
    
    \item The \bold{permutation-based test} is more expensive in terms of computational cost but it can be used whatever the sample size \eqn{n} is. The initial output samples (stored in the object of class \code{sensiHSIC}) are randomly permuted \eqn{B} times and the test statistic is recomputed as many times. This allows to simulate \eqn{B} observations of the test statistic under \eqn{H0} and to estimate the p-value in a non-parametric way. See Meynaoui (2019) for more details on how to correctly estimate the p-value in order to preserve the expected level of the test.
    
    \item The \bold{sequential permutation-based test} is a goal-oriented variant of the previous test. The main idea is to reduce the computational cost by stopping permutations as soon as the estimation of the p-value has sufficiently converged so that it can be compared to a reference threshold or be given a final ranking. See El Amri and Marrel (2022) for more details on how to implement this sequential approach for the three stopping criteria (namely \code{"ranking"}, \code{"screening"} or \code{"both"}).
    
    \item The \bold{Gamma test} is a parametric alternative to permutation-based tests when \eqn{n} is not large enough to resort to the asymptotic test. The permutation-based test reveals the test statistic under \eqn{H0} follows a unimodal distribution having significant positive skewness. Thus, it seems quite natural to estimate the p-value with a Gamma distribution, especially in view of the fact that the asymptotic distribution is properly approximated by this parametric family. See El Amri and Marrel (2021) for more details on how to estimate the parameters of the Gamma distribution with the method of moments. In particular, the first two moments of the test statistic under \eqn{H0} are computed thanks to the formulas that were initially provided in Kazi-Aoual et al. (1995).
    
  }

If a \bold{U-statistic} was used in \code{sensiHSIC}, the estimated value of \eqn{HSIC(Xi,Y)} may be negative. 

  \itemize{
    \item The \bold{asymptotic test} can no longer be conducted with a Gamma distribution (whose support is limited to \eqn{[0,+\infty[}). It is replaced by a Pearson III  distribution (which is a left-shifted Gamma distribution).
    \item The \bold{permutation-based test} and the \bold{sequential permutation-based test} can be applied directly.
    \item The \bold{Gamma test} has no longer any theoretical justification.
  }

}

\subsection{What about target and conditional \acronym{HSIC} indices?}{

In Marrel and Chabridon (2021), \acronym{HSIC} indices were adapted to \bold{target sensitivity analysis} (thus becoming \acronym{T-HSIC} indices) and to \bold{conditional sensitivity analysis} (thus becoming C-HSIC indices). Tests of independence can still be useful after estimating \acronym{T-HSIC} indices or \acronym{C-HSIC} indices.

  \itemize{
  
    \item For \acronym{T-HSIC} indices, the null hypothesis is \eqn{H0}: "\eqn{Xi} and \eqn{w(Y)} are independent" where \eqn{w} is the weight function selected in \code{target} and passed to the function \eqn{sensiHSIC}. Everything works just as for basic \acronym{HSIC} indices (apart from the fact that \eqn{w} is applied on the original output variable \eqn{Y}). Available test methods include \code{"Asymptotic"}, \code{"Permutation"}, \code{"Seq_Permutation"} and \code{"Gamma"} (for V-statistics only).
    
    \item For \acronym{C-HSIC} indices, the null hypothesis is \eqn{H0}: "\eqn{Xi} and \eqn{Y} are independent if the event described in \code{cond} occurs". In this specific context, testing conditional independence is only relevant if the weight function is an indicator function. For this reason, if conditional independence has to be tested, the user must select \code{type="indicTh"} in \code{cond} while calling the function \code{sensiHSIC}. Let us recall that only V-statistic estimators can be used for \acronym{C-HSIC} indices. As a result, available test methods include \code{"Asymptotic"}, \code{"Permutation"}, \code{"Seq_Permutation"} and \code{"Gamma"}.
    
  }
  
}

}









\references{

  El Amri, M. R. and Marrel, A. (2022), \emph{Optimized \acronym{HSIC}-based tests for sensitivity analysis: application to thermalhydraulic simulation of accidental scenario on nuclear reactor}, Quality and Reliability Engineering International, 38(3), 1386-1403.

  El Amri, M. R. and Marrel, A. (2021), \emph{More powerful \acronym{HSIC}-based independence tests, extension to space-filling designs and functional data}.
  \url{https://hal-cea.archives-ouvertes.fr/cea-03406956/}

  Fukumizu, K., Bach, F. R. and Jordan, M. I. (2004), \emph{Dimensionality reduction for supervised learning with reproducing kernel Hilbert spaces}, Journal of Machine Learning Research, 5(Jan), 73-99.
  
  Gretton, A., Fukumizu, K., Teo, C., Song, L., Scholkopf, B. and Smola, A. (2007), \emph{A kernel statistical test of independence}, Advances in Neural Information Processing Systems, 20.
  
  Kazi-Aoual, F., Hitier, S., Sabatier, R. and Lebreton, J. D. (1995), \emph{Refined approximations to permutation tests for multivariate inference}, Computational Statistics & Data Analysis, 20(6), 643-656.

  Marrel, A. and Chabridon, V. (2021), \emph{Statistical developments for target and conditional sensitivity analysis: application on safety studies for nuclear reactor}, Reliability Engineering & System Safety, 214, 107711.

  Meynaoui, A. (2019), \emph{New developments around dependence measures for sensitivity analysis: application to severe accident studies for generation IV reactors} (Doctoral dissertation, INSA de Toulouse).
  
  Sriperumbudur, B., Fukumizu, K. and Lanckriet, G. (2010), \emph{On the relation between universality, characteristic kernels and \acronym{RKHS} embedding of measures}, Proceedings of the 13th International Conference on Artificial Intelligence and Statistics (pp. 773-780). JMLR Workshop and Conference Proceedings.

}







\author{
  Sebastien Da Veiga, Amandine Marrel, Anouar Meynaoui, Reda El Amri and Gabriel Sarazin.
}










\seealso{
  \code{\link{sensiHSIC}, \link{weightTSA}}
}












\examples{
  \donttest{

# Test case: the Ishigami function.

n <- 20   # very few input-output samples
p <- 3    # nb of input variables

########################################
### PRELIMINARY SENSITIVITY ANALYSIS ###
########################################

X <- matrix(runif(n*p), n, p)
sensi <- sensiHSIC(model=ishigami.fun, X)
print(sensi)
plot(sensi)
title("GSA for the Ishigami function")

#############################
### TESTS OF INDEPENDENCE ###
#############################

test.asymp <- testHSIC(sensi)

test.perm <- testHSIC(sensi, test.method="Permutation")

test.seq.screening <- testHSIC(sensi, test.method="Seq_Permutation")

test.seq.ranking <- testHSIC(sensi, test.method="Seq_Permutation", 
                             seq.options=list(criterion="ranking"))

test.seq.both <- testHSIC(sensi, test.method="Seq_Permutation", 
                          seq.options=list(criterion="both"))

test.gamma <- testHSIC(sensi, test.method="Gamma")

# comparison of p-values

res <- rbind( t(as.matrix(test.asymp$pval)), t(as.matrix(test.perm$pval)), 
              t(as.matrix(test.seq.screening$pval)), t(as.matrix(test.seq.ranking$pval)),
              t(as.matrix(test.seq.both$pval)), t(as.matrix(test.gamma$pval)) )

rownames(res) <- c("asymp", "perm", "seq_perm_screening", 
                   "seq_perm_ranking", "seq_perm_both", "gamma")
res

# Conclusion: n is too small for the asymptotic test.
# Take n=200 and all four test methods will provide very close p-values.

#####################
### VISUALIZATION ###
#####################

# simulated values of HSIC indices under H0 (random permutations)
Hperm <- t(unname(test.perm$Hperm))

for(i in 1:p){
  
  # histogram of the test statistic under H0 (random permutations)
  
  title <- paste0("Histogram of S", i, " = HSIC(X", i, ",Y)")
  
  hist(Hperm[,i], probability=TRUE,
       nclass=70, main=title, xlab="", ylab="", col="cyan")
  
  # asymptotic Gamma distribution
  
  shape.asymp <- test.asymp$param[i, "shape"]
  scale.asymp <- test.asymp$param[i, "scale"]
  
  xx <- seq(0, max(Hperm[,i]), length.out=200)
  dens.asymp <- dgamma(xx, shape=shape.asymp, scale=scale.asymp)
  
  lines(xx, dens.asymp, lwd=2, col="darkorchid")
  
  # finite-sample Gamma distribution
  
  shape.perm <- test.gamma$param[i, "shape"]
  scale.perm <- test.gamma$param[i, "scale"]
  
  dens.perm <- dgamma(xx, shape=shape.perm, scale=scale.perm)
  
  lines(xx, dens.perm, lwd=2, col="blue")
  
  all.cap <- c("Asymptotic Gamma distribution", "Finite-sample Gamma distribution")
  all.col <- c("darkorchid", "blue")
  
  legend("topright", legend=all.cap, col=all.col, lwd=2, y.intersp=1.3)

}
 }
}
