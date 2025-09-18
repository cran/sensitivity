\name{testHSIC}
\alias{testHSIC}
\alias{print.testHSIC}
\alias{plot.testHSIC}




\title{Tests of Independence based on the Hilbert-Schmidt Independence Criterion (HSIC)}






\description{ 
\code{testHSIC} allows to test independence among all input-output pairs \eqn{(X_i, Y)} after a preliminary sensitivity analysis based on \acronym{HSIC} indices. \code{testHSIC} takes an object of class \code{sensiHSIC} (produced by a prior call to the function \code{sensiHSIC} that estimates \acronym{HSIC} indices) and it returns the estimated p-values after testing independence among all input-output pairs. For each input-output pair, having access to the p-value helps the user decide whether the null hypothesis \eqn{H_0}: "\eqn{X_i} and \eqn{Y} are independent" must be accepted or rejected. If all kernels selected in \code{sensiHSIC} are \bold{characteristic}, \eqn{H_0} can be rewritten "\eqn{HSIC(X_i, Y)=0}" and this paves the way to several test procedures.

Depending on the sample size, the selected \acronym{HSIC}-based sensitivity measure and the estimator type (either a \bold{U-statistic} or a \bold{V-statistic}), there are up to four different methods to test \eqn{H_0}. The \bold{asymptotic test} is recommended when the sample size \eqn{n} is around a few hundreds (or more). When \eqn{n} is smaller, a \bold{permutation-based test} must be considered instead. As a general rule, permutation-based tests can always be applied but a much heavier computational load is to be expected. However, if \acronym{HSIC} indices were initially estimated with V-statistics, the \bold{Gamma test} is a parametric method that offers an enticing tradeoff.

If the total-order \acronym{HSIC-ANOVA} indices were estimated during the preliminary sensitivity analysis, their numerators \eqn{T_i} also allow to test independence among all input-output pairs \eqn{(X_i, Y)}. For a given input variable \eqn{X_i}, let \eqn{T_i} denote the sum of \eqn{HSIC(X_i, Y)} and all \acronym{HSIC-ANOVA} interaction effects between \eqn{X_i} and other inputs. If all kernels selected in \code{sensiHSIC} are \bold{characteristic}, the null hypothesis \eqn{H_0}: "\eqn{X_i} and \eqn{Y} are independent" can be rewritten \eqn{H_0}: "\eqn{T_i=0}". The resulting test procedures prove to be more powerful, which means the probability of truly detecting dependence among \eqn{(X_i, Y)} is increased.

The test statistic associated to \eqn{T_i} can be either a \bold{U-statistic} or a \bold{V-statistic}. Several kinds of test procedures can then be considered. They are all inspired from what is done for \acronym{HSIC} indices. A \bold{permutation-based test} can always been applied but it becomes too expensive (from a computational viewpoint) when the sample size \eqn{n} is large. A much cheapier strategy consists in using the \bold{Gamma test} but it is limited to the specific case where \eqn{T_i} is estimated with a V-statistic.
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

  \item{sensi}{An object of class \code{"sensiHSIC"} which is produced by a prior call to the function \code{sensiHSIC}. In particular, \code{sensi} must contain objects named \code{"KX"} (3D-array filled with all input Gram matrices), \code{"KY"} (output Gram matrix), \code{"HSICXY"} (estimated \acronym{HSIC} indices) and \code{"estimator.type"} (either \code{"U-stat"} or \code{"V-stat"}). 
  
    \itemize{
      \item If the intended test procedure is based on total-order \acronym{HSIC-ANOVA} indices, \code{sensi} must also contain an object named \code{"TO.num"} (numerators of the estimated total-order \acronym{HSIC-ANOVA} indices).
      \item If the intended test procedure is based on \acronym{C-HSIC} indices, \code{sensi} must also contain objects named \code{"cond"} (list of options describing the conditioning event) and \code{"weights"} (conditioning weights after normalization).
    }
    
  }

  \item{test.method}{A string specifying the numerical procedure used to estimate the p-values. 
    
    \itemize{
      \item For a test of independence based on \acronym{HSIC} indices, available procedures include \code{"Asymptotic"} (asymptotic test), \code{"Permutation"} (permutation-based test), \code{"Seq_Permutation"} (sequential permutation-based test) and \code{"Gamma"} (Gamma test). 
      
        \itemize{
          \item If \code{sensi} contains V-statistics, the asymptotic test (resp. the Gamma test) is recommended for large (resp. small) sample sizes. Otherwise, permutation-based tests can be used as well.
          \item If \code{sensi} contains U-statistics, the Gamma test must not be employed. The asymptotic test is recommended for large sample sizes. Otherwise, permutation-based tests can be used as well. 
        }

      \item For a test of independence based on total-order \acronym{HSIC-ANOVA} indices, available procedures include \code{"Tot_Permutation"} (permutation-based test), \code{"Tot_Seq_Permutation"} (sequential permutation-based test) and \code{"Tot_Gamma"} (Gamma test). 
      
        \itemize{
          \item If \code{sensi} contains V-statistics, the Gamma test is recommended. Otherwise, permutation-based tests can be used as well.
          \item If \code{sensi} contains U-statistics, only permutation-based tests can be used.
        }

    }
    
  }
  
  \item{B}{Number of random permutations carried out on the output samples before the non-parametric estimation of p-values. Only relevant if \code{test.method} is \code{"Permutation"} or \code{"Tot_Permutation"}.}
  
  \item{seq.options}{A list of options guiding the sequential procedure. Only relevant if \code{test.method} is \code{"Seq_Permutation"} or \code{"Tot_Seq_Permutation"}.
  
    \itemize{
      \item \code{criterion} is a string specifying the stopping criterion. Available criteria include \code{"screening"} (permutations stop as soons as the estimated p-values have sufficiently converged so that they can be compared to the reference threshold \code{alpha}), \code{"ranking"} (permutations stop as soon as the estimated p-values have sufficiently converged so that they can be ranked) and \code{"both"} (permutations stop as soon as the two previous criteria are fulfilled).
      \item \code{alpha} is a scalar value (between \eqn{0} and \eqn{1}) specifying the type I error (probability of wrongly accepting \eqn{H_0}). Only relevant if \code{criterion} is \code{"screening"} or \code{"both"}.
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
  
  \item{family}{Only if \code{test.method} is \code{"Asymptotic"}, \code{"Gamma"} or \code{"Tot_Gamma"}. \cr
  A string indicating the parametric family used to estimate p-values.}
  
  \item{param}{Only if \code{test.method} is \code{"Asymptotic"}, \code{"Gamma"} or \code{"Tot_Gamma"}. \cr 
  A \eqn{2}-column (resp. \eqn{3}-column) matrix containing the parameters of the Gamma (resp. Pearson III) distributions used to estimate p-values.}
  
  \item{Hperm}{Only if \code{test.method} is \code{"Permutation"} or \code{"Tot_Permutation"}. \cr
  A \eqn{B}-column matrix containing simulated values of the test statistics after randomly permuting the output samples. Each column in \code{Hperm} corresponds to one random permutation.}
  
  \item{paths}{Only if \code{test.method} is \code{"Seq_Permutation"} or \code{"Tot_Seq_Permutation"}. \cr
  A matrix containing all estimated p-values over the sequential test procedure. The \eqn{i}-th row provides all estimates of the \eqn{i}-th p-value as the number of permutations increases. If one row ends with a sequence of missing values \code{NA}, it means permutations were stopped earlier for this input variable. This situation can only happen if \code{criterion="screening"} in the list \code{seq.options}.}

}










\details{

\subsection{Why and how to properly choose kernels?}{

For a given input-output pair of variables, the Hilbert-Schmidt independence criterion \acronym{(HSIC)} is a dissimilarity measure between the joint bivariate distribution and the product of marginal distributions. Dissimilarity between those two distributions is measured through the squared norm of the distance between their respective embeddings in a reproducing kernel Hilbert space (\acronym{RKHS}) that directly depends on the selected input kernel \eqn{K_{X_i}} and the selected output kernel \eqn{K_Y}.

It must always be kept in mind that this criterion allows to detect independence within the pair \eqn{(X_i, Y)} provided that the two kernels are \bold{characteristic}.  
  \itemize{
    \item If both kernels are characteristic, \eqn{H_0}: "\eqn{X_i} and \eqn{Y} are independent" is equivalent to \eqn{H_0}: "\eqn{HSIC(X_i, Y)=0}" and any estimator of \eqn{HSIC(X_i, Y)} emerges as a relevant test statistic.
    \item If they are not, testing \eqn{H_0}: "\eqn{HSIC(X_i, Y)=0}" is no longer sufficient for testing  \eqn{H_0}: "\eqn{X_i} and \eqn{Y} are independent".
  }
The reader is referred to Fukumizu et al. (2004) for the mathematical definition of a characteristic kernel and to Sriperumbur et al. (2010) for an overview of the major related results. 

Moreover, in order to achieve the \acronym{HSIC-ANOVA} decomposition, specific kernels called \bold{\acronym{ANOVA} kernels} must be used for the input variables. See Da Veiga (2021) for further details on \acronym{ANOVA} kernels, especially examples of ANOVA kernels and several methods to transform a basic kernel into an ANOVA kernel. Using \acronym{ANOVA} kernels is a necessary condition in order to define rigorously first-order and total-order \acronym{HSIC-ANOVA} indices.

In a nutshell, two situations must be distinguished:
  \enumerate{
    \item If the intended test of independence is built upon \eqn{HSIC(X_i, Y)}, all kernels must be \bold{characteristic}.
    \item If the intended test of independence is built upon \eqn{T_i}, all kernels must be \bold{characteristic} and \bold{\acronym{ANOVA} kernels} must be assigned to all input variables.
  }
  
Responsability for kernel selection is left to the user while calling the function \code{sensiHSIC}. Let us simply recall which kernels verify the above properties:

  \itemize{
  
    \item \bold{Characteristic kernel on any discrete set:} \code{"categ"} (categorical kernel).
  
    \item \bold{Characteristic kernels on \eqn{\mathbb{R}}:} \code{"laplace"} (exponential kernel), \code{"matern3"} (Matern \eqn{3/2} kernel), \code{"matern5"} (Matern \eqn{5/2} kernel) and \code{"rbf"} (Gaussian kernel). They remain characteristic when they are restricted to a compact domain \eqn{D} within \eqn{\mathbb{R}}.
    
    \item \bold{Characteristic kernels on \eqn{[0,1]}:} \code{"sobolev1"} (Sobolev kernel with smoothness parameter \eqn{r=1}) and \code{"sobolev2"} (Sobolev kernel with smoothness parameter \eqn{r=2}). In addition, the four abovementioned kernels (which are characteristic on \eqn{\mathbb{R}}) remain characteristic on \eqn{[0,1]}. Their \acronym{ANOVA} variants are also characteristic on \eqn{[0,1]}.
    
    \item \bold{Characteristic \acronym{ANOVA} kernel for any discrete distribution:} \code{"categ_anova"} (transformed categorical kernel).
    
    \item \bold{Characteristic \acronym{ANOVA} kernels for the uniform distribution on \eqn{[0,1]}:} \code{"sobolev1"} (Sobolev kernel with smoothness parameter \eqn{r=1}), \code{"sobolev2"} (Sobolev kernel with smoothness parameter \eqn{r=2}), \code{"laplace_anova"} (transformed exponential kernel), \code{"matern3_anova"} (transformed Matern \eqn{3/2} kernel), \code{"matern5_anova"} (transformed Matern \eqn{5/2} kernel) and \code{"rbf_anova"} (transformed Gaussian kernel).
    
  }
  
}

\subsection{Which test method is most appropriate?}{

First of all, the user must decide whether \eqn{H_0} has to be tested with \eqn{HSIC(X_i,Y)} or with \eqn{T_i}. It mainly depends on what is stored in the object \code{sensi} of class \code{"sensiHSIC"}. When both test statistics are available, it is strongly recommended to choose \eqn{T_i} because the associated test procedures proved to be more powerful in most cases. Unfortunately, in many situations (especially when the input variables \eqn{X_i} are not mutually independent), \eqn{T_i} cannot be estimated and \eqn{HSIC(X_i,Y)} is the only available option.

\bold{Tests of independence based on \acronym{HSIC} indices}

The test statistic for the pair \eqn{(X_i, Y)} is either the \bold{U-statistic} or the \bold{V-statistic} associated to \eqn{HSIC(X_i, Y)}.

If a \bold{V-statistic} was used in \code{sensiHSIC}, four different test methods can be considered.

  \itemize{
  
    \item The \bold{asymptotic test} can be used if the sample size \eqn{n} is large enough (at least a hundred of samples). The asymptotic distribution of the test statistic is approximated by a Gamma distribution whose parameters are estimated with the method of moments. See Gretton et al. (2007) for more details about how to estimate the first two moments of the asymptotic Gamma distribution.
    
    \item The \bold{permutation-based test} is more expensive in terms of computational cost but it can be used whatever the sample size \eqn{n} is. The initial output samples (stored in the object of class \code{sensiHSIC}) are randomly permuted \eqn{B} times and the test statistic is recomputed as many times. This allows to simulate \eqn{B} observations of the test statistic under \eqn{H_0} and to estimate the p-value in a non-parametric way. See Meynaoui (2019) for more details on how to correctly estimate the p-value in order to preserve the expected level of the test.
    
    \item The \bold{sequential permutation-based test} is a goal-oriented variant of the previous test. The main idea is to reduce the computational cost by stopping permutations as soon as the estimation of the p-value has sufficiently converged so that it can be compared to a reference threshold or be given a final ranking. See El Amri and Marrel (2022) for more details on how to implement this sequential approach for the three stopping criteria (namely \code{"ranking"}, \code{"screening"} or \code{"both"}).
    
    \item The \bold{Gamma test} is a parametric alternative to permutation-based tests when \eqn{n} is not large enough to resort to the asymptotic test. The permutation-based test reveals the test statistic under \eqn{H_0} follows a unimodal distribution having significant positive skewness. Thus, it seems quite natural to estimate the p-value with a Gamma distribution, especially in view of the fact that the asymptotic distribution is properly approximated by this parametric family. See El Amri and Marrel (2021) for more details on how to estimate the parameters of the Gamma distribution with the method of moments. In particular, the first two moments under \eqn{H_0} of the test statistic are computed thanks to the formulas that were initially provided in Kazi-Aoual et al. (1995).
    
  }

If a \bold{U-statistic} was used in \code{sensiHSIC}, the estimated value of \eqn{HSIC(X_i,Y)} may be negative. 

  \itemize{
    \item The \bold{asymptotic test} can no longer be conducted with a Gamma distribution (whose support is limited to \eqn{[0,+\infty[}). It is replaced by a Pearson III  distribution (which is a left-shifted Gamma distribution).
    \item The \bold{permutation-based tests} can be applied directly.
    \item The \bold{Gamma test} has no longer any theoretical justification.
  }

\bold{Tests of independence based on total-order \acronym{HSIC-ANOVA} indices}

The test statistic for the pair \eqn{(X_i, Y)} is either the \bold{U-statistic} or the \bold{V-statistic} associated to \eqn{T_i}.

If a \bold{V-statistic} was used in \code{sensiHSIC}, three different test methods can be considered.

  \itemize{
  
    \item The \bold{permutation-based test} is a simple adaptation of what is done for \eqn{HSIC(X_i,Y)}. A slightly different permutation scheme is repeated \eqn{B} times and the test statistic is recomputed as many times. This allows to simulate \eqn{B} observations of the test statistic under \eqn{H_0} and to estimate the p-value in a non-parametric way. This test may be very expensive if \eqn{n} is large.
    
    \item The \bold{sequential permutation-based test} is the goal-oriented variant of the previous test. It is a simple adaptation of what is done for \eqn{HSIC(X_i,Y)}. It allows to reduce as much as possible the computational load induced while recomputing the test statistic for a series of random permutations.
    
    \item The \bold{Gamma test} is a parametric procedure that may be used when the computational load of permutation-based tests becomes too heavy (which happens when the sample size \eqn{n} is large). It is an extension of what is done for \eqn{HSIC(X_i,Y)}. In order to compute the first two moments of the test statistic under \eqn{H_0}, the analytical formulas proposed in Kazi-Aoual et al. (1995) were altered so as to comply with less accomodating assumptions.
    
  }

If a \bold{U-statistic} was used in \code{sensiHSIC}, the estimated value of \eqn{T_i} may be negative. Only \bold{permutation-based tests} can be considered.

}

\subsection{What about target and conditional \acronym{HSIC} indices?}{

In Marrel and Chabridon (2021), \acronym{HSIC} indices were adapted to \bold{target sensitivity analysis} (thus becoming \acronym{T-HSIC} indices) and to \bold{conditional sensitivity analysis} (thus becoming C-HSIC indices). Tests of independence can still be useful after estimating \acronym{T-HSIC} indices or \acronym{C-HSIC} indices.

  \itemize{
  
    \item For \acronym{T-HSIC} indices, the null hypothesis is \eqn{H_0}: "\eqn{X_i} and \eqn{w(Y)} are independent" where \eqn{w} is the weight function selected in \code{target} and passed to the function \code{sensiHSIC}. Everything works just as for basic \acronym{HSIC} indices (apart from the fact that \eqn{w} is applied on the original output variable \eqn{Y}). Available test methods include \code{"Asymptotic"}, \code{"Permutation"}, \code{"Seq_Permutation"} and \code{"Gamma"} (for V-statistics only).
    
    \item For \acronym{C-HSIC} indices, the null hypothesis is \eqn{H_0}: "\eqn{X_i} and \eqn{Y} are independent if the event described in \code{cond} occurs". In this specific context, testing conditional independence is only relevant if the weight function is an indicator function. For this reason, if conditional independence has to be tested, the user must select \code{type="indicTh"} in \code{cond} while calling the function \code{sensiHSIC}. Let us recall that only V-statistic estimators can be used for \acronym{C-HSIC} indices. As a result, available test methods include \code{"Asymptotic"}, \code{"Permutation"}, \code{"Seq_Permutation"} and \code{"Gamma"}.
    
  }
  
\bold{Target sensitivity analysis} can also be considered with total-order \acronym{HSIC-ANOVA} indices (when the context allows to define them, which requires among other things mutually independent input variables). If their numerators were estimated during the prior call to \code{sensiHSIC}, the null hypothesis \eqn{H_0}: "\eqn{X_i} and \eqn{w(Y)} are independent" can be tested with \eqn{T_i}. Associated test methods include \code{"Tot_Permutation"}, \code{"Tot_Seq_Permutation"} and \code{"Tot_Gamma"} (for V-statistics only). 

\bold{Conditional sensitivity analysis} cannot be considered with total-order \acronym{HSIC-ANOVA} indices because the conditioning operation is feared to introduce statistical dependence among input variables, which forbids using \acronym{HSIC-ANOVA} indices.
  
}

}










\references{

  Da Veiga, S. (2021), \emph{Kernel-based ANOVA decomposition and Shapley effects: application to global sensitivity analysis}, arXiv preprint arXiv:2101.05487.

  El Amri, M. R. and Marrel, A. (2022), \emph{Optimized \acronym{HSIC}-based tests for sensitivity analysis: application to thermalhydraulic simulation of accidental scenario on nuclear reactor}, Quality and Reliability Engineering International, 38(3), 1386-1403.

  El Amri, M. R. and Marrel, A. (2021), \emph{More powerful \acronym{HSIC}-based independence tests, extension to space-filling designs and functional data}.

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
# Test case: the Matyas function.

n <- 300  # nb of samples
p <- 3    # nb of input variables (including 1 dummy variable)

########################################
### PRELIMINARY SENSITIVITY ANALYSIS ###
########################################

X <- matrix(runif(n*p), n, p)
sensi <- sensiHSIC(model=matyas.fun, X, 
                   kernelX="sobolev1", anova=list(obj="both", is.uniform=TRUE))
print(sensi)
plot(sensi)
title("GSA for the Matyas function")

#############################
### TESTS OF INDEPENDENCE ###
#############################

### HSIC indices ###

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

### total-order HSIC-ANOVA indices ###

test.tot.perm <- testHSIC(sensi, test.method="Tot_Permutation")

test.tot.seq.screening <- testHSIC(sensi, test.method="Tot_Seq_Permutation")

test.tot.seq.ranking <- testHSIC(sensi, test.method="Tot_Seq_Permutation", 
                             seq.options=list(criterion="ranking"))

test.tot.seq.both <- testHSIC(sensi, test.method="Tot_Seq_Permutation", 
                          seq.options=list(criterion="both"))

test.tot.gamma <- testHSIC(sensi, test.method="Tot_Gamma")

res <- rbind( t(as.matrix(test.tot.perm$pval)), 
              t(as.matrix(test.tot.seq.screening$pval)), 
              t(as.matrix(test.tot.seq.ranking$pval)),
              t(as.matrix(test.tot.seq.both$pval)), 
              t(as.matrix(test.tot.gamma$pval)) )

rownames(res) <- c("tot_perm", "tot_seq_perm_screening", 
                   "tot_seq_perm_ranking", "tot_seq_perm_both", "tot_gamma")
res

#####################
### VISUALIZATION ###
#####################

# simulated values of HSIC indices under H0 (random permutations)
Hperm <- t(unname(test.perm$Hperm))

# simulated values of total-order HSIC-ANOVA indices under H0 (random permutations)
tot.Hperm <- t(unname(test.tot.perm$Hperm))

for(i in 1:p){
  
  ttl <- paste0("First-order and total-order HSIC-ANOVA indices for X", i)
  
  ######################################
  ### FIRST-ORDER HSIC-ANOVA INDICES ###
  ######################################
  
  # histogram of the test statistic under H0 (random permutations)
  
  hist(Hperm[,i], probability=TRUE,
       nclass=70, main=ttl, xlab="", ylab="", col="cyan")
  
  xx <- seq(0, max(tot.Hperm[,i]), length.out=200)
  
  # asymptotic Gamma distribution
  
  shape.asymp <- test.asymp$param[i, "shape"]
  scale.asymp <- test.asymp$param[i, "scale"]
  
  dens.asymp <- dgamma(xx, shape=shape.asymp, scale=scale.asymp)
  
  lines(xx, dens.asymp, lwd=2, col="darkorchid")
  
  # finite-sample Gamma distribution
  
  shape.perm <- test.gamma$param[i, "shape"]
  scale.perm <- test.gamma$param[i, "scale"]
  
  dens.perm <- dgamma(xx, shape=shape.perm, scale=scale.perm)
  
  lines(xx, dens.perm, lwd=2, col="blue")
  
  ######################################
  ### TOTAL-ORDER HSIC-ANOVA INDICES ###
  ######################################
  
  # histogram of the test statistic under H0 (random permutations)
  
  hist(tot.Hperm[,i], probability=TRUE, add=TRUE,
       nclass=70, xlab="", ylab="", col="deeppink")
  
  # finite-sample Gamma distribution
  
  shape.tot.perm <- test.tot.gamma$param[i, "shape"]
  scale.tot.perm <- test.tot.gamma$param[i, "scale"]
  
  dens.tot.perm <- dgamma(xx, shape=shape.tot.perm, scale=scale.tot.perm)
  
  lines(xx, dens.tot.perm, lwd=2, col="red")
  
  ### legend ###
  
  txt.1 <- paste0("Histogram of HSIC(X", i, ",Y)")
  txt.11 <- "Asymptotic Gamma distribution"
  txt.12 <- "Finite-sample Gamma distribution"
  
  txt.2 <- paste0("Histogram of T", i, " = HSIC(X,Y) - HSIC(X", 
                  paste(setdiff(1:p, i), collapse=""), ",Y)")
  txt.21 <- "Finite-sample Gamma distribution"
  
  all.cap <- c(txt.1, txt.11, txt.12, txt.2, txt.21)
  all.col <- c("cyan", "darkorchid", "blue", "deeppink", "red")
  all.lwd <- c(7, 2, 2, 7, 2)
  
  legend("topright", legend=all.cap, col=all.col, lwd=all.lwd, y.intersp=1.3)

}
 }
}
