\name{sensiHSIC}
\alias{sensiHSIC}
\alias{tell.sensiHSIC}
\alias{print.sensiHSIC}
\alias{plot.sensiHSIC}














\title{Sensitivity Indices based on the Hilbert-Schmidt Independence Criterion (HSIC)}
















\description{ 

\code{sensiHSIC} allows to conduct \bold{global sensitivity analysis (\acronym{GSA})} in many different contexts thanks to several sensitivity measures based on the \bold{Hilbert-Schmidt independence criterion (\acronym{HSIC})}. The so-called \acronym{HSIC} sensitivity indices depend on the kernels which are affected to the input variables \eqn{Xi} as well as on the kernel which is affected to the output object \eqn{Y}. For each random entity, a reproducing kernel Hilbert space (\acronym{RKHS}) is associated to the chosen kernel and allows to represent probability distributions in an appropriate function space. The influence of \eqn{Xi} on \eqn{Y} is then measured through the distance between the joint probability distribution (true impact of \eqn{Xi} on \eqn{Y} in the numerical model) and the product of marginal distributions (no impact of \eqn{Xi} on \eqn{Y}) after embedding those distributions into a bivariate \acronym{RKHS}. Such a GSA approach has three main advantages:
  \itemize{
    \item The input variables \eqn{Xi} may be correlated.
    \item Any kind of mathematical object is supported (provided that a kernel function is available).
    \item Accurate estimation is possible even in presence of very few data (no more than a hundred of input-output samples).
  }
In \code{sensiHSIC}, each input variable \eqn{Xi} is expected to be scalar (either discrete or continous). On the contrary, a much wider collection of mathematical objects are supported for the output variable \eqn{Y}. In particular, \eqn{Y} may be:

  \itemize{
  
    \item A \bold{scalar output} (either discrete or continous). If so, one single kernel family is selected among the kernel collection.
    
    \item A \bold{low-dimensional vector output}. If so, a kernel is selected for each output variable and the final output kernel is built by \bold{tensorization}.
    
    \item A \bold{high-dimensional vector output} or a \bold{functional output}. In this case, the output data must be seen as time series observations. Three different methods are proposed.
    
    \enumerate{
      \item A preliminary \bold{dimension reduction} may be performed. In order to achieve this, a \bold{principal component analysis (\acronym{PCA})} based on the empirical covariance matrix helps identify the first terms of the Kharunen-Loeve expansion. The final output kernel is then built in the reduced subspace where the functional data are projected.
      \item The \bold{dynamic time warping (\acronym{DTW})} algorithm may be combined with a translation-invariant kernel. The resulting \acronym{DTW}-based output kernel is well-adapted to measure similarity between two given time series.
      \item The \bold{global alignment kernel (\acronym{GAK})} may be directly applied on the functional data. Unlike the \acronym{DTW} kernel, it was specifically designed to deal with time series and to measure the impact of one input variable on the shape of the output curve.
    }
    
  }

Many variants of the original \acronym{HSIC} indices are now available in \code{sensiHSIC}.
  
  \itemize{

    \item \bold{Normalized \acronym{HSIC} indices (\acronym{R2-HSIC})} \cr
    The original \acronym{HSIC} indices defined in Gretton et al. (2005) may be hard to interpret because they do not admit a universal upper bound. A first step to overcome this difficulty was enabled by Da Veiga (2015) with the definition of the \acronym{R2-HSIC} indices. The resulting sensitivity indices can no longer be greater than \eqn{1}.
  
    \item \bold{Target \acronym{HSIC} indices (\acronym{T-HSIC})} \cr
    They were thought by Marrel and Chabridon (2021) to meet the needs of \bold{target sensitivity analysis (\acronym{TSA})}. The idea is to measure the impact of each input variable \eqn{Xi} on a specific part of the output distribution (for example the upper tail). To achieve this, a weight function \eqn{w} is applied on \eqn{Y} before computing \acronym{HSIC} indices.
  
    \item \bold{Conditional \acronym{HSIC} indices (\acronym{C-HSIC})} \cr
    They were thought by Marrel and Chabridon (2021) to meet the needs of \bold{conditional sensitivity analysis (\acronym{CSA})}. The idea is to measure the impact of each input variable \eqn{Xi} on \eqn{Y} when a specific event occurs. This conditioning event is detected on the output variable \eqn{Y} and its amplitude is taken into account thanks to a weight function \eqn{w}.
  
    \item \bold{\acronym{HSIC-ANOVA} indices} \cr
    To improve the interpretability of \acronym{HSIC} indices, Da Veiga (2021) came up with an \bold{ANOVA-like decomposition} that allows to establish a strict separation of main effects and interaction effects in the \acronym{HSIC} paradigm. The first-order and total-order \acronym{HSIC-ANOVA} indices are then defined just as the first-order and total-order Sobol' indices. However, this framework only holds if two major assumptions are verified: the input variables \eqn{Xi} must be mutually independent and all input kernels must belong to the very restrained class of \acronym{ANOVA} kernels. 
    
  }

As most sensitivity measures, \acronym{HSIC} indices allow to rank the input variables \eqn{Xi} according to the influence they have on the output variable \eqn{Y}. They can also be used for a screening purpose, that is to distinguish between truly influential input variables and non-influential input variables. The user who is interested in this topic is invited to consult the documentation of the function \code{testHSIC}. 

}













\usage{
sensiHSIC(model = NULL, X, target = NULL, cond = NULL, 
          kernelX = "rbf", paramX = NA,
          kernelY = "rbf", paramY = NA,
          estimator.type = "V-stat",
          nboot = 0, conf = 0.95,
          anova = list(obj = "no", is.uniform = TRUE),
          sensi = NULL, 
          save.GM = list(KX = TRUE, KY = TRUE), \dots)
          
\method{tell}{sensiHSIC}(x, y = NULL, \dots)

\method{print}{sensiHSIC}(x, \dots)

\method{plot}{sensiHSIC}(x, ylim = c(0, 1), \dots)
}
















\arguments{

  \item{model}{A function, or a statistical model with a \code{predict} method. It defines the input-output model that needs to be studied.}
  
  \item{X}{A \eqn{n}-by-\eqn{p} matrix containing all input samples. It comprises \eqn{n} joint observations of the \eqn{p} input variables.
  
    \itemize{
      \item If the user is only wanting to estimate \acronym{HSIC} indices or \acronym{R2-HSIC} indices, the input variables can be correlated.
      \item If the user is also wanting to estimate \acronym{HSIC-ANOVA} indices, the input variables have to be mutually independent.
    }

  }
  
  \item{target}{A list of options to perform TSA. At least, \code{target} must contain an option named \code{"c"}. For other options, there exist default assignments.
    \itemize{
      \item \code{type} is a string specifying the weight function. Available choices include \code{"indicTh"}, \code{"zeroTh"}, \code{"logistic"} and \code{"exp1side"}. Default value is \code{"exp1side"}.
        \itemize{
          \item \code{"indicTh"} and \code{"zeroTh"} only depend on a threshold parameter. 
          \item \code{"logistic"} and \code{"exp1side"} depend on both a threshold parameter and a smoothness parameter.
        }
      \item \code{c} is a scalar value specifying the threshold parameter.
      \item \code{upper} is a boolean indicating whether the target region is located above (\code{TRUE}) or below (\code{FALSE}) the threshold parameter \code{c}. Only relevant when \code{type} is \code{"indicTh"}, \code{"zeroTh"} or \code{"exp1side"}. Default value is \code{TRUE}. 
      \item \code{param} is a scalar value specifying the smoothness parameter. Only relevant when \code{type} is \code{"logistic"} or \code{"exp1side"}. Default value is \code{1}.
    }
  }
  
  \item{cond}{A list of options to perform CSA. At least, \code{cond} must contain an option named \code{"c"}. For other options, there exist default assignments.
    \itemize{
      \item \code{type} is a string specifying the weight function. Available choices include \code{"indicTh"}, \code{"zeroTh"}, \code{"logistic"} and \code{"exp1side"}. Default value is \code{"exp1side"}.
        \itemize{
          \item \code{"indicTh"} and \code{"zeroTh"} only depend on a threshold parameter. 
          \item \code{"logistic"} and \code{"exp1side"} depend on both a threshold parameter and a smoothness parameter.
        }
      \item \code{c} is a scalar value specifying the threshold parameter.
      \item \code{upper} is a boolean indicating whether the conditioning region is located above (\code{TRUE}) or below (\code{FALSE}) the threshold parameter \code{c}. Only relevant when \code{type} is \code{"indicTh"}, \code{"zeroTh"} or \code{"exp1side"}. Default value is \code{TRUE}. 
      \item \code{param} is a scalar value specifying the smoothness parameter. Only relevant if \code{type} is \code{"logistic"} or \code{"exp1side"}. Default value is \code{1}.
    }
  }
  
  \item{kernelX}{A string or a vector of \eqn{p} strings that specifies how to choose input kernels.
  
    \itemize{
      \item If only one string is provided, the associated kernel is affected to all inputs.
      \item For dimension-wise kernel selection, a vector of \eqn{p} strings must be provided.
    }
    
  For each input variable, available choices include \code{"categ"} (categorical kernel), \code{"dcov"} (covariance kernel of the fractional Brownian motion), \code{"invmultiquad"} (inverse multiquadratic kernel), \code{"laplace"} (exponential kernel), \code{"linear"} (dot-product kernel), \code{"matern3"} (Matern \eqn{3/2} kernel), \code{"matern5"} (Matern \eqn{5/2} kernel), \code{"raquad"} (rationale quadratic kernel), \code{"rbf"} (Gaussian kernel), \code{"sobolev1"} (Sobolev kernel with smoothness parameter \eqn{r=1}) and \code{"sobolev2"} (Sobolev kernel with smoothness parameter \eqn{r=2}).

  In addition, let us assume that all input variables are uniformly distributed on \eqn{[0,1]}. Under this assumption, the kernels \code{"laplace"}, \code{"matern3"}, \code{"matern5"} and \code{"rbf"} can be easily transformed into \acronym{ANOVA} kernels. The resulting kernels are respectively called \code{"laplace_anova"}, \code{"matern3_anova"}, \code{"matern5_anova"} and \code{"rbf_anova"}.
  
    \itemize{
      \item One-parameter kernels: \code{"categ"}, \code{"dcov"}, \code{"invmultiquad"}, \code{"laplace"}, \code{"laplace_anova"}, \code{"matern3"}, \code{"matern3_anova"}, \code{"matern5"}, \code{"matern5_anova"}, \code{"raquad"}, \code{"rbf"} and \code{"rbf_anova"}.
      \item Parameter-free kernels: \code{"linear"}, \code{"sobolev1"} and \code{"sobolev2"}. 
    }
  
  }
    
  \item{paramX}{A scalar value or a vector of \eqn{p} values with input kernel parameters.
  
    \itemize{
      \item If \code{paramX=NA}, input kernel parameters are computed automatically with rules of thumb.
      \item If \code{paramX} is a scalar value, it is affected to all input kernels.
      \item For dimension-wise kernel parametrization, a vector of \eqn{p} values must be provided. If \code{kernelX} combines one-parameter kernels and parameter-free kernels, \code{NA} must be specified for parameter-free kernels.
    }
    
  } 
  
  \item{kernelY}{A string, a vector of \eqn{q} strings or a list of options that specifies how to construct the output kernel. Regardless of its mathematical nature, the model output must be envisioned as a \eqn{q}-dimensional random vector.
  
  To deal with a \bold{scalar output} or a \bold{low-dimensional vector output}, it is advised to select one kernel per output dimension and to tensorize all selected kernels. In this case, \code{kernelY} must be a string or a vector of \eqn{q} strings.

  \itemize{
    \item If only one string is provided, the associated kernel is repeated \eqn{q} times.
    \item For dimension-wise kernel selection, a vector of \eqn{q} strings must be provided.
  }
  
  Have a look at \code{kernelX} for an exhaustive list of available kernels.
  
  To deal with a \bold{high-dimensional vector output} or a \bold{functional output}, it is advised to reduce dimension or to use a dedicated kernel. In this case, \code{kernelY} must be specified as a list of options. At least, \code{kernelY} must contain an option named \code{"method"}. For other options, there exist default assignments.
  
  \itemize{
    \item \code{method} is a string indicating the strategy used to construct the output kernel. Available choices include \code{"PCA"} (dimension reduction through principal component analysis), \code{"DTW"} (dynamic type warping) and \code{"GAK"} (global alignment kernel).
  }
  
    \enumerate{
  
      \item If \code{method="PCA"}, the following options may also be specified:
  
        \itemize{
          \item \code{data.centering} is a boolean indicating whether the input samples must be centered before performing the preliminary PCA. Default value is \code{TRUE}.
          \item \code{data.scaling} is a boolean indicating whether the input samples must be scaled before performing the preliminary PCA. Default value is \code{TRUE}.
          \item \code{fam} is a string specifying the input kernel which is applied on principal components. Available choices only include \code{"dcov"}, \code{"invmultiquad"}, \code{"laplace"}, \code{"linear"}, \code{"matern3"}, \code{"matern5"}, \code{"raquad"} and \code{"rbf"}. Default value is \code{"rbf"}.
          \item \code{expl.var} is a scalar value (between \eqn{0} and \eqn{1}) specifying the expected percentage of output variance that must be explained by \acronym{PCA}. Default value is \code{0.95}.
          \item \code{PC} is the expected number of principal components in \acronym{PCA}. Default value is \code{NA}.
          \item \code{combi} is a string indicating how the final output kernel is built in the reduced subspace. Available options include \code{"sum"} or \code{"prod"}. The chosen kernel in \code{fam} is applied on all principal components before summation (if \code{"sum"}) or tensorization (if \code{"prod"}).
          \item \code{position} is a string indicating whether weights have to be involved in the construction of the final output kernel in the reduced subspace. Available choices include \code{"nowhere"} (no weights), \code{"intern"} (weights applied on principal components) or \code{"extern"} (weights applied on kernels). Default value is \code{"intern"}.
        }
        
      \bold{Remark:} \code{expl.var} and \code{PC} are conflicting options. Only one of both needs to be specified and the other one must be set to \code{NA}. If both are specified, \code{expl.var} is prioritized. If both are set to \code{NA}, \code{expl.var} is then set to its default value.
  
      \item If \code{method="DTW"}, the following option may also be specified:
  
        \itemize{
          \item \code{fam} is a string specifying the translation-invariant kernel which is combined with \acronym{DTW}. Available choices only include \code{"invmultiquad"}, \code{"laplace"}, \code{"matern3"}, \code{"matern5"}, \code{"raquad"} and \code{"rbf"}. Default value is \code{"rbf"}.
        }
        
      \item If \code{method="GAK"}, there is no other option to specify.
  
    }

  }
  
  \item{paramY}{A scalar value or a vector of values with output kernel parameters.
    
    \itemize{
      \item If \code{paramY=NA}, output kernel parameters are computed automatically with rules of thumb.
    }
    
    In other cases, \code{paramY} must be specified in agreement with \code{kernelY}. 
    
    \bold{Case 1:} \code{kernelY} is a string or a vector of \eqn{q} strings.
    
    \code{paramY} must be a scalar value or a vector of \eqn{q} values with output kernel parameters.
    
      \itemize{
        \item If \code{paramY} is a scalar value, it is affected to all output kernels.
        \item For dimension-wise kernel parametrization, a vector of \eqn{q} values must be provided. If \code{kernelY} combines one parameter kernels and parameter-free kernels, \code{NA} must be specified for parameter-free kernels.
      }
    
    \bold{Case 2:} \code{kernelY} is a list of options with \code{method="PCA"}.
    
    \code{paramY} must be set to \code{NA} because the parameters involved in the final output kernel are computed automatically. Their optimal tuning depends on the reduced subspace given by \acronym{PCA}.
    
    \bold{Case 3:} \code{kernelY} is a list of options with \code{method="DTW"}.
    
    \code{paramY} must be set to \code{NA}.

    \bold{Case 4:} \code{kernelY} is a list of options with \code{method="GAK"}.
    
    \code{paramY} must be a vector of \eqn{2} values. If the user only wants to specify one parameter, the other one must be set to \code{NA}. The two parameters correspond to the arguments \code{sigma} and \code{window.size} in the function \code{gak} from the package \code{dtwclust}. However, automatical computation (specified by \code{paramY=NA}) is strongly advised for this kind of output kernel.

  }
  
  \item{estimator.type}{A string specifying the kind of estimator used for HSIC indices. Available choices include \code{"U-stat"} (U-stastics) and \code{"V-stat"} (V-statistics). U-statistics are unbiased estimators. V-statistics are biased estimators but they become unbiased asymptotically. In the specific case of \acronym{HSIC} indices, V-statistics are non-negative estimators and they offer more flexibility for further test procedures (see \code{testHSIC}). Both kinds of estimators can be computed with complexity \eqn{O(n^2)} where \eqn{n} denotes the sample size.}

  \item{nboot}{Number of bootstrap replicates.}
  
  \item{conf}{A scalar value (between \eqn{0} and \eqn{1}) specifying the level of confidence intervals.}

  \item{anova}{A list of parameters to achieve an \acronym{ANOVA}-like decomposition based on HSIC indices. At least, \code{anova} must contain an option named \code{"obj"}. For other options, there exist default assignments.
  
    \itemize{
    
      \item \code{obj} is a string specifying which kinds of \acronym{HSIC-ANOVA} indices are expected. Available choices include \code{"no"} (\code{anova} is disabled), \code{"FO"} (first-order only), \code{"TO"} (total-order only) and \code{"both"} (first-order and total-order). 
      \item \code{is.uniform} is a boolean indicating whether the samples stored in \code{X} come from random variables that are uniformly distributed on \eqn{[0,1]}. Let us recall that \acronym{HSIC-ANOVA} indices can only be computed by means of \acronym{ANOVA} kernels. Among available kernels, only \code{"laplace_anova"}, \code{"matern3_anova"}, \code{"matern5_anova"}, \code{"rbf_anova"}, \code{"sobolev1"} and \code{"sobole2"} verify this constraint (provided that all input variables are uniformly distributed on \eqn{[0,1]}).
      
      \itemize{
        \item If \code{is.uniform=TRUE}, it is checked that each input data stored in \eqn{X} actually lies in \eqn{[0,1]}. If this condition is not verified, an error is returned.
        \item If \code{is.uniform=FALSE}, non-parametric rescaling (based on empirical distribution functions) is operated. 
      }
      
    }
    
  }
  
  \item{sensi}{An object of class \code{"sensiHSIC"} resulting from a prior call to the hereby function. If an object of class \code{"sensiHSIC"} is indeed provided, the following happens:
  
    \itemize{
      \item If \code{sensi} contains an object named \code{"KX"}, it is extracted from \code{sensi} and the input Gram matrices (required to estimate \acronym{HSIC} indices) are not built from \code{X}, \code{kernelX} and \code{paramX}.
      \item If \code{sensi} contains an object named \code{"KY"}, it is extracted from \code{sensi} and the output Gram matrix (required to estimate \acronym{HSIC} indices) is not built from \code{model}, \code{kernelY} and \code{paramY}. 
    }
    
  }

  \item{save.GM}{A list of parameters indicating whether Gram matrices have to be saved. The list \code{save.GM} must contain options named \code{"KX"} and \code{"KY"}.
    
    \itemize{
      
      \item \code{KX} is a boolean indicating whether the input Gram matrices have to be saved.
      \item \code{KY} is a boolean indicating whether the output Gram matrix has to be saved.
      
    }
    
  }

  \item{x}{An object of class \code{"sensiHSIC"} storing the state of the sensitivity study (parameters, data, estimates).}
  
  \item{y}{A \eqn{n}-by-\eqn{q} matrix containing all output samples. It comprises \eqn{n} observations of the \eqn{q} output variables.}
  
  \item{ylim}{A vector of two values specifying the \eqn{y}-coordinate plotting limits.}

  \item{\dots}{Any other arguments for \code{model} which are passed unchanged each time \code{model} is called.}
      
}














\value{

  \code{sensiHSIC} returns a list of class \code{"sensiHSIC"}. It contains all the input arguments detailed before, except \code{sensi} which is not kept. It must be noted that some of them might have been altered, corrected or completed.
  
  \item{kernelX}{A vector of \eqn{p} strings with input kernels.}
  
  \item{paramX}{A vector of \eqn{p} values with input kernel parameters. For each one-parameter kernel, a real number is returned. It is either the original value (if correct), a corrected value (if not) or the default value (computed from a rule of thumb when \code{NA} is specified). For each parameter-free kernel, \code{NA} is returned.}
  
  \item{kernelY}{A vector of \eqn{q} strings or a list of options that specifies how the output kernel was constructed. In the case where \code{kernelY} is a list of options with \code{method="PCA"}, \code{kernelY} contains additional information resulting from \acronym{PCA}.
    \itemize{
      \item If \code{kernelY} initally contained an option named \code{"expl.var"}, \code{kernelY} now also contains an option named \code{"PC"} that provides the associated number of principal components.
      \item If \code{kernelY} initially contained an option named \code{"PC"}, \code{kernelY} now also contains an option named \code{"expl.var"} that provides the associated percentage of output variance that is explained by \acronym{PCA}.
      \item If \code{kernelY} initally contained an option named \code{"position"} that was set to \code{"intern"} or \code{"extern"}, \code{kernelY} now contains an option named \code{"ratios"} that provides the weights used to combine kernels in the reduced subspace given by \acronym{PCA}.
    }
  }
  
  \item{paramY}{A vector of values with output kernel parameters.
  
    \bold{Case 1:} \code{kernelY} is a list of \eqn{q} strings.
    
    \code{paramY} is a vector of \code{q} values. For each one-parameter kernel, a real number is returned. It is either the original value (if correct), a corrected value or the default value (computed with a rule of thumb if \code{NA} was initially specified). For each parameter-free kernel, \code{NA} is returned.
  
    \bold{Case 2:} \code{kernelY} is a list of options with \code{method="PCA"}.
    
    \code{paramY} is a vector of \code{PC} values. For this method, let us recall that all kernels belong to the same family which is specified by an option named \code{"fam"} within \code{kernelY}. For each dimension in the reduced subspace, the kernel parameter is computed (with a rule of thumb) from the corresponding principal component. If the kernel in \code{fam} is parameter-free, \code{paramY} is a vector where \code{NA} is repeated \code{PC} times.
    
    \bold{Case 3:} \code{kernelY} is a list of options with \code{method="DTW"}.
    
    \code{paramY} remains equal to \code{NA}.
    
    \bold{Case 4:} \code{kernelY} is a list of options with \code{method="GAK"}.
    
    \code{paramY} is a vector of \eqn{2} values. For each parameter, the returned value is either the original value (if correct), a corrected value or the default value (computed with a rule of thumb if \code{NA} was initially specified). 
  
  }
  
  More importantly, the list of class \code{"sensiHSIC"} contains all expected results (output samples, sensitivity measures and conditioning weights).
  
  \item{call}{The matched call.}
  
  \item{y}{A \eqn{n}-row matrix containing all output samples. The \eqn{i}-th row in \code{y} is obtained from the \eqn{i}-th row in \code{X} after computing the model response. If \code{target} is passed to \code{sensiHSIC}, output samples in \code{y} are obtained after applying consecutively \code{model} and the specified weight function.}
  
  \item{HSICXY}{The estimated \acronym{HSIC} indices.}
  
  \item{S}{The estimated \acronym{R2-HSIC} indices (also called normalized \acronym{HSIC} indices).}
  
  \item{weights}{Only if \code{cond} is passed to \code{sensiHSIC}. \cr
  A vector of \eqn{n} values containing all conditioning weights. In the \acronym{CSA} context, the conditioning factor is defined by \eqn{w(Y)/E[w(Y)]}. See Marrel and Chabridon (2021) for further explanations.}
  
  Depending on what is specified in \code{anova}, the list of class \code{"sensiHSIC"} may also contain the following objects:

  \item{FO}{The estimated first-order \acronym{HSIC-ANOVA} indices.}
  
  \item{TO}{The estimated total-order \acronym{HSIC-ANOVA} indices.}
  
  \item{TO.num}{The estimated numerators of total-order \acronym{HSIC-ANOVA} indices.}
  
  \item{denom}{The estimated common denominator of \acronym{HSIC-ANOVA} indices.}
  
}












\details{

Let \eqn{(Xi,Y)} be an input-output pair. The kernels assigned to \eqn{Xi} and \eqn{Y} are respectively denoted by \eqn{Ki} and \eqn{KY}. 

For many global sensitivity measures, the influence of \eqn{Xi} on \eqn{Y} is measured in the light of the probabilistic dependence that exists within the input-output pair \eqn{(Xi,Y)}. For this, a dissimilarity measure is applied between the joint probability distribution (true impact of \eqn{Xi} and \eqn{Y} in the numerical model) and the product of marginal distributions (no impact of \eqn{Xi} on \eqn{Y}). For instance, Borgonovo's sensitivity measure is built upon the total variation distance between those two probability distributions. See Borgonovo and Plischke (2016) for further details.

The \acronym{HSIC}-based sensitivity measure can be understood in this way since the index \eqn{HSIC(Xi,Y)} results from the application of the \bold{Hilbert-Schmidt independence criterion (\acronym{HSIC})} on the pair \eqn{(Xi,Y)}. This criterion is nothing but a special kind of dissimilarity measure between the joint probability distribution and the product of marginal distributions. This dissimilarity measure is called the \bold{maximum mean discrepancy (MMD)} and its definition relies on the selected kernels \eqn{Ki} and \eqn{KY}. According to the theory of reproducing kernels, every kernel \eqn{K} is related to a \bold{reproducing kernel Hilbert space (\acronym{RKHS})}.Then, if \eqn{K} is affected to a random variable \eqn{Z}, any probability distribution describing the random behavior of \eqn{Z} may be represented within the induced \acronym{RKHS}. In this setup, the dissimilarity between the joint probability distribution and the product of marginal distributions is then measured through the squared norm of their images into the bivariate \acronym{RKHS}. The user is referred to Gretton et al. (2006) for additional details on the mathematical construction of \acronym{HSIC} indices.

In practice, it may be difficult to understand how \eqn{HSIC(Xi,Y)} measures dependence within \eqn{(Xi,Y)}. An alternative definition relies on the concept of \bold{feature map}. Let us recall that the value taken by a kernel function can always be seen as the scalar product of two \bold{feature functions} lying in a \bold{feature space}. Definition 1 in Gretton et al. (2005) introduces \eqn{HSIC(Xi,Y)} as the Hilbert-Schmidt norm of a covariance-like operator between random features. For this reason, having access to the input and output feature maps may help identify the dependence patterns captured by \eqn{HSIC(Xi,Y)}.

Kernels must be chosen very carefully. There exists a wide variety of kernels but only a few f them meet the needs of \acronym{GSA}. As \eqn{HSIC(Xi,Y)} is supposed to be a dependence measure, it must be equal to \eqn{0} if and only if \eqn{Xi} and \eqn{Y} are independent. A sufficient condition to enable this equivalence is to take two characteristic kernels. The reader is referred to Fukumizu et al. (2004) for the mathematical definition of a characteristic kernel and to Sriperumbur et al. (2010) for an overview of the major related results. In particular:
  \itemize{
    \item The Gaussian kernel, the Laplace kernel, the Matern \eqn{3/2} kernel and the Matern \eqn{5/2} kernel (all defined on \eqn{R^2}) are \bold{characteristic}.
    \item The transformed versions of the four abovementioned kernels (all defined on \eqn{[0,1]^2}) are \bold{characteristic}.
    \item All Sobolev kernels (defined on \eqn{[0,1]^2}) are \bold{characteristic}.
    \item The categorical kernel (defined on any discrete probability space) is \bold{characteristic}.
  }
  
Lemma 1 in Gretton et al. (2005) provides a third way of defining \eqn{HSIC(Xi,Y)}. Since the associated formula is only based on three expectation terms, the corresponding estimation procedures are very simple and they do not ask for a large amount of input-output samples to be accurate. Two kinds of estimators may be used for \eqn{HSIC(Xi,Y)}: the \bold{V-statistic estimator} (which is non negative, biased and asymptotically unbiased) or the \bold{U-statistic estimator} (unbiased). For both estimators,  the computational complexity is \eqn{O(n^2)} where \eqn{n} is the sample size.

The user must always keep in mind the key steps leading to the estimation of \eqn{HSIC(Xi,Y)}:
  \itemize{
    \item Input samples are simulated and the corresponding output samples are computed with the numerical model.
    \item An input kernel \eqn{Ki} and an output kernel \eqn{KY} are selected.
    \item \bold{In case of target sensitivity analysis:} output samples are transformed by means of a weight function \eqn{w}.
    \item The input and output Gram matrices are constructed.
    \item \bold{In case of conditional sensitivity analysis:} conditioning weights are computed by means of a weight function \eqn{w}.
    \item The final estimate is computed. It depends on the selected estimator type (either a U-statistic or a V-statistic).
  }
  
\subsection{Kernel functions for random variables}{

All what follows is written for a scalar output \eqn{Y} but the same is true for any scalar input \eqn{Xi}. 

Let \eqn{D} denote the support of the output probability distribution. A kernel is a symmetric and positive definite function defined on the domain \eqn{D}. Different kernel families are available in \code{sensiHSIC}. 
  \itemize{
    \item To deal with continuous probability distributions on \eqn{R}, one can use:
    \itemize{
      \item The covariance kernel of the fractional Browian motion (\code{"dcov"}), the inverse multiquadratic kernel (\code{"invmultiquad"}), the exponential kernel (\code{"laplace"}), the dot-product kernel (\code{"linear"}), the Matern \eqn{3/2} kernel (\code{"matern3"}), the Matern \eqn{5/2} kernel (\code{"matern5"}), the rationale quadratic kernel (\code{"raquad"}) and the Gaussian kernel (\code{"rbf"}).
    }
    \item To deal with continuous probability distributions on \eqn{[0,1]}, one can use:
    \itemize{
      \item Any of the abovementioned kernel (restricted to \eqn{[0,1]}).
      \item The transformed exponential kernel (\code{"laplace_anova"}), the transformed Matern \eqn{3/2} kernel (\code{"matern3_anova"}), the transformed Matern \eqn{5/2} kernel (\code{"matern5_anova"}), the transformed Gaussian kernel (\code{"rbf_anova"}), the Sobolev kernel with smoothness parameter \eqn{r=1} (\code{"sobolev1"}) and the Sobolev kernel with smoothness parameter \eqn{r=2} (\code{"sobolev2"}).
    }
    \item To deal with any discrete probability distribution, the categorical kernel (\code{"categ"}) must be used.
  }
  
Two kinds of kernels must be distinguished:

  \itemize{
  
    \item \bold{Parameter-free kernels:} the dot-product kernel (\code{"linear"}), the Sobolev kernel with smoothness parameter \eqn{r=1} (\code{"sobolev1"}) and the Sobolev kernel with smoothness parameter \eqn{r=2} (\code{"sobolev2"}).
    
    \item \bold{One-parameter kernels:} the categorical kernel (\code{"categ"}), the covariance kernel of the fractional Brownian motion kernel (\code{"dcov"}), the inverse multiquadratic kernel (\code{"invmultiquad"}), the exponential kernel (\code{"laplace"}), the transformed exponential kernel (\code{"laplace_anova"}), the Matern \eqn{3/2} kernel (\code{"matern3"}), the transformed Matern \eqn{3/2} kernel (\code{"matern3_anova"}), the Matern \eqn{5/2} kernel (\code{"matern5"}), the transformed Matern \eqn{5/2} kernel (\code{"matern5_anova"}), the rationale quadratic kernel (\code{"raquad"}), the Gaussian kernel (\code{"rbf"}) and the transformed Gaussian kernel (\code{"rbf_anova"}).
    
  }

A major issue related to one-parameter kernels is how to set the parameter. It mainly depends on the role played by the parameter in the kernel expression.

  \itemize{
  
    \item For translation-invariant kernels and their \acronym{ANOVA} variants (that is all one-parameter kernels except \code{"categ"} and \code{"dcov"}), the parameter may be interpreted as a correlation length (or a scale parameter). The rule of thumb is to compute the empirical standard deviation of the provided samples.
    
    \item For the covariance kernel of the fractional Brownian motion (\code{"dcov"}), the parameter is an exponent. Default value is \eqn{1}.
    
    \item For the categorical kernel (\code{"categ"}), the parameter has no physical sense. It is just a kind of binary encoding.
      \itemize{
        \item \eqn{0} means the user wants to use the basic categorical kernel.
        \item \eqn{1} means the user wants to use the weighted variant of the categorical kernel.
      }
      
  }

}

\subsection{How to deal with a low-dimensional vector output?}{

Let us assume that the output vector \eqn{Y} is composed of \eqn{q} random variables \eqn{Y1,...,Yq}.

A kernel \eqn{Kj} is affected to each output variable \eqn{Yj} and this leads to embed the \eqn{j}-th output probability distribution in a \acronym{RKHS} denoted by \eqn{Hj}. Then, the \bold{tensorization} of \eqn{H1,...,Hq} allows to build the final \acronym{RKHS}, that is the \acronym{RKHS} where the \eqn{q}-variate output probability distribution describing the overall random behavior of \eqn{Y} will be embedded. In this situation:
  \itemize{
    \item The final output kernel is the tensor product of all output kernels.
    \item The final output Gram matrix is the Hadamard product of all output Gram matrices.
  }
Once the final output Gram matrix is built, \acronym{HSIC} indices can be estimated, just as in the case of a scalar output.

}

\subsection{How to deal with a high-dimensional vector output or a functional output?}{

In \code{sensiHSIC}, three different methods are proposed in order to compute \acronym{HSIC}-based sensitivity indices in presence of functional outputs.

\bold{Dimension reduction}

This approach was initially proposed by Da Veiga (2015). The key idea is to approximate the random functional output by the first terms of its \bold{Kharunen-Loeve expansion}. This can be achived with a \bold{principal component analysis (PCA)} that is carried out on the empirical covariance matrix.
  \itemize{
    \item The eigenvectors (or \bold{principal directions}) allow to approximate the (deterministic) functional terms involved in the Kharunen-Loeve decomposition.
    \item The eigenvalues allow to determine how many principal directions are sufficient in order to accurately represent the random function by means of its truncated Kharunen-Loeve expansion. The key idea behind dimension reduction is to keep as few principal directions as possible while preserving a prescribed level of explained variance.
  }
The \bold{principal components} are the coordinates of the functional output in the low-dimensional subspace resulting from \acronym{PCA}. There are computed for all output samples (time series observations). See Le Maitre and Knio (2010) for more detailed explanations.

The last step consists in constructing a kernel in the reduced subspace. One single kernel family is selected and affected to all principal directions. Moreover, all kernel parameters are computed automatically (with appropriate rules of thumb). Then, several strategies may be considered.

  \itemize{
  
    \item The initial method described in Da Veiga (2015) is based on a direct tensorization. One can also decide to sum kernels.
    
    \item This approach was improved by El Amri and Marrel (2021). For each principal direction, a weight coefficient (equal the ratio between the eigenvalue and the sum of all selected eigenvalues) is computed.
      \itemize{
        \item The principal components are multiplied by their respective weight coefficients before summing kernels or tensorizing kernels.
        \item The kernels can also be directly applied on the principal components before being linearly combined according to the weight coefficients.
      }

  }
In sensiHSIC, all these strategies correspond to the following specifications in \code{kernelY}:
  \itemize{
    
    \item \bold{Direct tensorization:}
    \code{kernelY=list(method="PCA", combi="prod", position="nowhere")}
    
    \item \bold{Direct sum:}
    \code{kernelY=list(method="PCA", combi="sum", position="nowhere")}
    
    \item \bold{Rescaled tensorization:}
    \code{kernelY=list(method="PCA", combi="prod", position="intern")}
    
    \item \bold{Rescaled sum:} 
    \code{kernelY=list(method="PCA", combi="sum", position="intern")}
    
    \item \bold{Weighted linear combination:}
    \code{kernelY=list(method="PCA", combi="sum", position="extern")}
    
  }
  
\bold{Dynamic Time Warping (\acronym{DTW})}

The \acronym{DTW} algorithm developed by Sakoe and Chiba (1978) can be combined with a translation-invariant kernel in order to create a kernel function for times series. The resulting \acronym{DTW}-based output kernel is well-adapted to measure similarity between two given time series.

Suitable translation-invariant kernels include the inverse multiquadratic kernel (\code{"invmultiquad"}), the exponential kernel (\code{"laplace"}), the Matern \eqn{3/2} kernel (\code{"matern3"}), the Matern \eqn{5/2} kernel (\code{"matern5"}), the rationale quadratic kernel (\code{"raquad"}) and the Gaussian kernel (\code{"rbf"}).

The user is warned against the fact that \acronym{DTW}-based kernels are not positive definite functions. As a consequence, many theoretical properties do not hold anymore for \acronym{HSIC} indices.

For faster computations, \code{sensiHSIC} is using the function \code{dtw_dismat} from the package \code{incDTW}.

\bold{Global Alignment Kernel (\acronym{GAK})}

Unlike \acronym{DTW}-based kernels, the \acronym{GAK} is a positive definite function. This time-series kernel was originally introduced in Cuturi et al. (2007) and further investigated in Cuturi (2011). It was used to compute \acronym{HSIC} indices on a simplified compartmental epidemiological model in Da Veiga (2021).

For faster computations, \code{sensiHSIC} is using the function \code{gak} from the package \code{dtwclust}. 

In \code{sensiHSIC}, two \acronym{GAK}-related parameters may be tuned by the user with \code{paramY}. They exactly correspond to the arguments \code{sigma} and \code{window.size} in the function \code{gak}.

}

\subsection{About normalized \acronym{HSIC} indices (\acronym{R2-HSIC})}{

No doubt interpretability is the major drawback of \acronym{HSIC} indices. This shortcoming led Da Veiga (2021) to introduce a normalized version of \eqn{HSIC(Xi,Y)}. The so-called R2-HSIC index is thus defined as the ratio between \eqn{HSIC(Xi,Y)} and the square root of a normalizing constant equal to \eqn{HSIC(Xi,Xi)*HSIC(Y,Y)}. 

This normalized sensitivity measure is inspired from the \bold{distance correlation measure} proposed by Szekely et al. (2007) and the resulting sensitivity indices are easier to interpret since they all fall in the interval \eqn{[0,1]}.

}

\subsection{About target \acronym{HSIC} indices (\acronym{T-HSIC})}{

T-HSIC indices were designed by Marrel and Chabridon (2021) for \acronym{TSA}. They are only defined for a scalar output. Vector and functional outputs are not supported. The main idea of \acronym{TSA} is to measure the influence of each input variable \eqn{Xi} on a modified version of \eqn{Y}. To do so, a preliminary mathematical transform \eqn{w} (called the \bold{weight function}) is applied on \eqn{Y}. The collection of \acronym{HSIC} indices is then estimated with respect to \eqn{w(Y)}. Here are two examples of situations where \acronym{TSA} is particularly relevant:
  \itemize{
    \item How to measure the impact of \eqn{Xi} on the upper values taken by \eqn{Y} (for example the values above a given threshold \eqn{T})?
      \itemize{
        \item To answer this question, one may take \eqn{w(Y)=Y*1_{Y>T}} (\bold{zero-thresholding}). \cr
        This can be specified in \code{sensiHSIC} with \code{target=list(c=T, type="zeroTh", upper=TRUE)}.
      }
    \item How to measure the influence of \eqn{Xi} on the occurrence of the event \eqn{{Y>T}}?
      \itemize{
        \item To answer this question, one may take \eqn{w(Y)=1_{Y<T}} (\bold{indicator-thresholding}). \cr
        This can be specified in \code{sensiHSIC} with \code{target=list(c=T, type="indicTh", upper=FALSE)}.
      }
  }

In Marrel and Chabridon (2021), the two situations described above are referred to as \bold{"hard thresholding"}. To avoid using discontinuous weight functions, \bold{"smooth thresholding"} may be used instead.
  \itemize{
    \item Spagnol et al. (2019): logistic transformation on both sides of the threshold \eqn{T}.
    \item Marrel and Chabridon (2021): exponential transformation above or below the threshold \eqn{T}.
  }
These two smooth relaxation functions depend on a tuning parameter that helps control smoothness. For further details, the user is invited to consult the documentation of the function \code{weightTSA}.

\bold{Remarks:}

\itemize{

  \item When \code{type="indicTh"} (\bold{indicator-thesholding}), \eqn{w(Y)} becomes a binary random variable. Accordingly, the output kernel selected in \code{kernelY} must be the categorical kernel.
  
  \item In the spirit of \acronym{R2-HSIC} indices, \acronym{T-HSIC} indices can be normalized. The associated normalizing constant is equal to the square root of \eqn{HSIC(Xi,Xi)*HSIC(w(Y),w(Y))}.
  
  \item \acronym{T-HSIC} indices can be very naturally combined with the \acronym{HSIC-ANOVA} decomposition proposed by Da Veiga (2021). As a consequence, the arguments \code{target} and \code{anova} in \code{sensiHSIC} can be enabled simultaneously. Compared with basic \acronym{HSIC} indices, there are three main differences: the input variables must be mutually independent, \acronym{ANOVA} kernels must be used for all input variables and the output of interest is \eqn{w(Y)}.
  
  \item \acronym{T-HSIC} indices can be very naturally combined with the tests of independence proposed in \code{testHSIC}. In this context, the null hypothesis is \eqn{H0}: "\eqn{Xi} and \eqn{w(Y)} are independent".
  
}

}

\subsection{About conditional \acronym{HSIC} indices (\acronym{C-HSIC})}{

\acronym{C-HSIC} indices were designed by Marrel and Chabridon (2021) for \acronym{CSA}. They are only defined for a scalar output. Vector and functional outputs are not supported. The idea is to measure the impact of each input variable \eqn{Xi} on \eqn{Y} when a specific event occurs. This conditioning event is defined on \eqn{Y} thanks to a \bold{weight function} \eqn{w}. In order to compute the conditioning weights, \eqn{w} is applied on the output samples and an empirical normalization is carried out (so that the overall sum of conditioning weights is equal to \eqn{1}). The conditioning weights are then combined with the simulated Gram matrices in order to estimate \acronym{C-HSIC} indices. All formulas can be found in Marrel and Chabridon (2021). Here is an exemple of a situation where \acronym{CSA} is particularly relevant:
  \itemize{
    \item  Let us imagine that the event \eqn{{Y>T}} coincides with a system failure. \cr
    How to measure the influence of \eqn{Xi} on \eqn{Y} when failure occurs?
    \itemize{
      \item To answer this question, one may take \eqn{w(Y) = 1_{Y>T}} (\bold{indicator-thresholding}). \cr
      This can be specified in \code{sensiHSIC} with \code{cond=list(c=T, type="indicTh", upper=TRUE)}.
    }
  }
The three other weight functions proposed for TSA (namely \code{"zeroTh"}, \code{"logistic"} and \code{"exp1side"}) can also be used but the role they play is less intuitive to understand. See Marrel and Chabridon (2021) for better explanations.

\bold{Remarks:}

  \itemize{
  
    \item Unlike what is pointed out for \acronym{TSA}, when \code{type="thresholding"}, the output of interest \eqn{Y} remains a continuous random variable. The categorical kernel is thus inappropriate. A continuous kernel must be used instead.
    
    \item In the spirit of \acronym{R2-HSIC} indices, \acronym{C-HSIC} indices can be normalized. The associated normalizing constant is equal to the square root of \eqn{C-HSIC(Xi,Xi)*C-HSIC(Y,Y)}.

    \item Only V-statistics are supported to estimate C-HSIC indices. The reason is because the normalized version of C-HSIC indices cannot always be estimated with U-statistics. In particular, the estimates of \eqn{C-HSIC(Xi,Xi)*C-HSIC(Y,Y)} may be negative.
    
    \item \acronym{C-HSIC} indices cannot be combined with the \acronym{HSIC-ANOVA} decomposition proposed in Da Veiga (2021). In fact, the conditioning operation is feared to introduce statistical dependence among input variables, which forbids using \acronym{HSIC-ANOVA} indices. As a consequence, the arguments \code{cond} and \code{anova} in \code{sensiHSIC} cannot be enabled simultaneously.
    
    \item \acronym{C-HSIC} indices can harly be combined with the tests of inpendence proposed in \code{testHSIC}. This is only possible if \code{type="indicTh"}. In this context, the null hypothesis is \eqn{H0}: "\eqn{Xi} and \eqn{Y} are independent if the event described in \code{cond} occurs".
    
  }

}

\subsection{About \acronym{HSIC-ANOVA} indices}{

In comparison with \acronym{HSIC} indices, \acronym{R2-HSIC} indices are easier to interpret. However, in terms of interpretability, Sobol' indices remain much more convenient since they can be understood as shares of the total output variance. Such an interpretation is made possible by the Hoeffding decomposition, also known as \acronym{ANOVA} decomposition.

It was proved in Da Veiga (2021) that an \acronym{ANOVA}-like decomposition can be achived for \acronym{HSIC} indices under certain conditions:
  \itemize{
    \item The input variables must be mutually independent (which was not required to compute all other kinds of \acronym{HSIC} indices).
    \item \bold{\acronym{ANOVA} kernels} must be assigned to all input variables.
  }
This \acronym{ANOVA} setup allows to establish a strict separation between main effects and interaction effects in the \acronym{HSIC} sense. The first-order and total-order \acronym{HSIC-ANOVA} indices are then defined in the same fashion than first-order and total-order Sobol' indices. It is worth noting that the \acronym{HSIC-ANOVA} normalizing constant is equal to \eqn{HSIC(X,Y)} and is thus different from the one used for \acronym{R2-HSIC} indices.

For a given probability measure \eqn{P}, an \acronym{ANOVA} kernel \eqn{K} is a kernel that can rewritten \eqn{1+k} where \eqn{k} is an orthogonal kernel with respect to \eqn{P}. Among the well-known parametric families of probability distributions and kernel functions, there are very few examples of orthogonal kernels. One example is given by \bold{Sobolev kernels} when there are matched with the uniform probability measure on [0,1]. See Wahba et al. (1995) for further details on Sobolev kernels.

Moreover, several strategies to construct orthogonal kernels from non-orthogonal kernels are recalled in Da Veiga (2021). One of them consists in translating the feature map so that the resulting kernel becomes centered at the prescribed probability measure \eqn{P}. This can be done analytically for some basic kernels (Gaussian, exponential, Matern \eqn{3/2} and Matern \eqn{5/2}) when \eqn{P} is the uniform measure on \eqn{[0,1]}. See Section 9 in Ginsbourger et al. (2016) for the corresponding formulas.

In \code{sensiHSIC}, \acronym{ANOVA} kernels are only available for the uniform probability measure on \eqn{[0,1]}. This includes the Sobolev kernel with parameter \eqn{r=1} (\code{"sobolev1"}), the Sobolev kernel with parameter \eqn{r=2} (\code{"sobolev2"}), the transformed Gaussian kernel (\code{"rbf_anova"}), the transformed exponential kernel (\code{"laplace_anova"}), the transformed Matern \eqn{3/2} kernel (\code{"matern3_anova"}) and the transformed Matern \eqn{5/2} kernel (\code{"matern5_anova"}).

As explained above, the \acronym{HSIC-ANOVA} indices can only be computed if all input variables are uniformly distributed on \eqn{[0,1]}. Because of this limitation, a preliminary reformulation is needed if the \acronym{GSA} problem includes other kinds of input probability distributions. The \bold{probability integral transform (PIT)} must be applied on each input variable \eqn{Xi}. In addition, all quantile functions must be encapsulated in the numerical model, which may lead to reconsider the way \code{model} is specified. In \code{sensiHSIC}, if \code{check=TRUE} is selected in \code{anova}, it is checked that all input samples lie in \eqn{[0,1]}. If this is not the case, a non-parametric rescaling (based on empirical distribution functions) is operated.

\acronym{HSIC-ANOVA} indices can be used for \acronym{TSA}. The only difference with \acronym{GSA} is the use of a weight function \eqn{w}. On the contrary, CSA cannot be conducted with \acronym{HSIC-ANOVA} indices. Indeed, the conditioning operation is feared to introduce statistical independence among the input variables, which prevents using the \acronym{HSIC-ANOVA} approach.

}

}

\references{

  Borgonovo, E. and Plischke, E. (2016), \emph{Sensitivity analysis: a review of recent advances}, European Journal of Operational Research, 248(3), 869-887.

  Cuturi, M., Vert, J. P., Birkenes, O. and Matsui, T. (2007), \emph{A kernel for time series based on global alignments}, 2007 IEEE International Conference on Acoustics, Speech and Signal Processing-ICASSP'07 (Vol. 2, pp. II-413), IEEE.

  Cuturi, M. (2011), \emph{Fast global alignment kernels}, Proceedings of the 28th International Conference on Machine Learning (ICML-11) (pp. 929-936).

  Da Veiga, S. (2015), \emph{Global sensitivity analysis with dependence measures}, Journal of Statistical Computation and Simulation, 85(7), 1283-1305.
  
  Da Veiga, S. (2021). \emph{Kernel-based \acronym{ANOVA} decomposition and Shapley effects: application to global sensitivity analysis}, arXiv preprint arXiv:2101.05487.
  
  El Amri, M. R. and Marrel, A. (2021), \emph{More powerful \acronym{HSIC}-based independence tests, extension to space-filling designs and functional data}.
  \url{https://hal-cea.archives-ouvertes.fr/cea-03406956/}

  Fukumizu, K., Bach, F. R. and Jordan, M. I. (2004), \emph{Dimensionality reduction for supervised learning with reproducing kernel Hilbert spaces}, Journal of Machine Learning Research, 5(Jan), 73-99.
  
  Ginsbourger, D., Roustant, O., Schuhmacher, D., Durrande, N. and Lenz, N. (2016), \emph{On \acronym{ANOVA} decompositions of kernels and Gaussian random field paths}, Monte Carlo and Quasi-Monte Carlo Methods (pp. 315-330), Springer, Cham.
  
  Gretton, A., Bousquet, O., Smola, A., and Scholkopf, B. (2005), \emph{Measuring statistical dependence with Hilbert-Schmidt norms}, International Conference on Algorithmic Learning Theory (pp. 63-77), Springer, Berlin, Heidelberg.
  
  Gretton, A., Borgwardt, K., Rasch, M., Scholkopf, B. and Smola, A. (2006), \emph{A kernel method for the two-sample-problem}, Advances in Neural Information Processing Systems, 19.

  Le Maitre, O. and Knio, O. M. (2010), \emph{Spectral methods for uncertainty quantification with applications to computational fluid dynamics}, Springer Science & Business Media.
  
  Marrel, A. and Chabridon, V. (2021), \emph{Statistical developments for target and conditional sensitivity analysis: application on safety studies for nuclear reactor}, Reliability Engineering & System Safety, 214, 107711.
  
  Sakoe, H. and Chiba, S. (1978), \emph{Dynamic programming algorithm optimization for spoken word recognition}, IEEE International Conference on Acoustics, Speech and Signal, 26(1), 43-49.
  
  Spagnol, A., Riche, R. L. and Veiga, S. D. (2019), \emph{Global sensitivity analysis for optimization with variable selection}, SIAM/ASA Journal on Uncertainty Quantification, 7(2), 417-443.
  
  Sriperumbudur, B., Fukumizu, K. and Lanckriet, G. (2010), \emph{On the relation between universality, characteristic kernels and \acronym{RKHS} embedding of measures}, Proceedings of the 13th International Conference on Artificial Intelligence and Statistics (pp. 773-780). JMLR Workshop and Conference Proceedings.
  
  Szekely, G. J., Rizzo, M. L. and Bakirov, N. K. (2007), \emph{Measuring and testing dependence by correlation of distances}, The Anals of Statistics, 35(6), 2769-2794.
  
  Wahba, G., Wang, Y., Gu, C., Klein, R. and Klein, B. (1995), \emph{Smoothing spline \acronym{ANOVA} for exponential families, with application to the Wisconsin Epidemiological Study of Diabetic Retinopathy: the 1994 Neyman Memorial Lecture}, The Annals of Statistics, 23(6), 1865-1895.

}









\author{
  Sebastien Da Veiga, Amandine Marrel, Anouar Meynaoui, Reda El Amri and Gabriel Sarazin.
}










\seealso{
\code{\link{testHSIC}, \link{weightTSA}}
}











\examples{
 \donttest{

############################
### HSIC indices for GSA ###
############################

# Test case 1: the Friedman function
# --> 5 input variables

### GSA with a given model ###

n <- 800
p <- 5
X <- matrix(runif(n*p), n, p)

kernelX <- c("rbf", "rbf", "laplace", "laplace", "sobolev1")
paramX <- c(0.2, 0.3, 0.4, NA, NA)

# kernel for X1: Gaussian kernel with given parameter 0.2
# kernel for X2: Gaussian kernel with given parameter 0.3
# kernel for X3: exponential kernel with given parameter 0.4
# kernel for X4: exponential kernel with automatic computation of the parameter
# kernel for X5: Sobolev kernel (r=1) with no parameter

kernelY <- "raquad"
paramY <- NA 

sensi <- sensiHSIC(model=friedman.fun, X,
                   kernelX=kernelX, paramX=paramX, 
                   kernelY=kernelY, paramY=paramY)

print(sensi)
plot(sensi)
title("GSA for the Friedman function")

### GSA with given data ###

Y <- friedman.fun(X)
sensi <- sensiHSIC(model=NULL, X,
                   kernelX=kernelX, paramX=paramX, 
                   kernelY=kernelY, paramY=paramY)
tell(sensi, y=Y)

print(sensi)

### GSA from a prior object of class "sensiHSIC" ###

new.sensi <- sensiHSIC(model=friedman.fun, X,
                       kernelX=kernelX, paramX=paramX, 
                       kernelY=kernelY, paramY=paramY,
                       estimator.type="U-stat", 
                       sensi=sensi,
                       save.GM=list(KX=FALSE, KY=FALSE))

print(new.sensi)

# U-statistics are computed without rebuilding all Gram matrices.
# Those Gram matrices are not saved a second time.

##################################
### HSIC-ANOVA indices for GSA ###
##################################

# Test case 2: the Matyas function with Gaussian input variables
# --> 3 input variables (including 1 dummy variable)

n <- 10^3
p <- 2

X <- matrix(rnorm(n*p), n, p)

# The Sobolev kernel (with r=1) is used to achieve the HSIC-ANOVA decomposition.
# Both first-order and total-order HSIC-ANOVA indices are expected.

### AUTOMATIC RESCALING ###

kernelX <- "sobolev1"
anova <- list(obj="both", is.uniform=FALSE)

sensi.A <- sensiHSIC(model=matyas.fun, X, kernelX=kernelX, anova=anova)

print(sensi.A)
plot(sensi.A)
title("GSA for the Matyas function")

### PROBLEM REFORMULATION ###

U <- matrix(runif(n*p), n, p)
new.matyas.fun <- function(U){ matyas.fun(qnorm(U)) }

kernelX <- "sobolev1"
anova <- list(obj="both", is.uniform=TRUE)

sensi.B <- sensiHSIC(model=new.matyas.fun, U, kernelX=kernelX, anova=anova)

print(sensi.B)

####################################
### T-HSIC indices for target SA ###
####################################

# Test case 3: the Sobol function
# --> 8 input variables

n <- 10^3
p <- 8

X <- matrix(runif(n*p), n, p)

kernelY <- "categ"
target <- list(c=0.4, type="indicTh")

sensi <- sensiHSIC(model=sobol.fun, X, kernelY=kernelY, target=target)

print(sensi)
plot(sensi)
title("TSA for the Sobol function")

#########################################
### C-HSIC indices for conditional SA ###
#########################################

# Test case 3: the Sobol function
# --> 8 input variables

n <- 10^3
p <- 8

X <- matrix(runif(n*p), n, p)

cond <- list(c=0.2, type="exp1side", upper=FALSE)

sensi <- sensiHSIC(model=sobol.fun, X, cond=cond)

print(sensi)
plot(sensi)
title("CSA for the Sobol function")

##########################################
### How to deal with discrete outputs? ###
##########################################

# Test case 4: classification of the Ishigami output
# --> 3 input variables
# --> 3 categories

classif <- function(X){
  
  Ytemp <- ishigami.fun(X) 
  Y <- rep(NA, n)
  Y[Ytemp<0] <- 0
  Y[Ytemp>=0 & Ytemp<10] <- 1                
  Y[Ytemp>=10] <- 2  
  
  return(Y)
  
}

###

n <- 10^3
p <- 3

X <- matrix(runif(n*p, -pi, pi), n, p)

kernelY <- "categ"
paramY <- 0

sensi <- sensiHSIC(model=classif, X, kernelY=kernelY, paramY=paramY)
print(sensi)
plot(sensi)
title("GSA for the classified Ishigami function")

############################################
### How to deal with functional outputs? ###
############################################

# Test case 5: the arctangent temporal function
# --> 3 input variables (including 1 dummy variable)

n <- 500
p <- 3

X <- matrix(runif(n*p,-7,7), n, p)

### with a preliminary dimension reduction by PCA ###

kernelY <- list(method="PCA", 
                data.centering=TRUE, data.scaling=TRUE,
                fam="rbf", expl.var=0.95, combi="sum", position="extern")

sensi <- sensiHSIC(model=atantemp.fun, X, kernelY=kernelY)

print(sensi)
plot(sensi)
title("PCA-based GSA for the arctangent temporal function")

### with a kernel based on dynamic time warping ###

kernelY <- list(method="DTW", fam="rbf")

sensi <- sensiHSIC(model=atantemp.fun, X, kernelY=kernelY)

print(sensi)
plot(sensi)
title("DTW-based GSA for the arctangent temporal function")


\donttest{
### with the global alignment kernel ###

kernelY <- list(method="GAK")

sensi <- sensiHSIC(model=atantemp.fun, X, kernelY=kernelY)

print(sensi)
plot(sensi)
title("GAK-based GSA for the arctangent temporal function")
}
  }
}
