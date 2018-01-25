\name{multiselect}
\alias{multiselect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selecting Combinations of Predictors by Leveraging Multiple AUCs for an Ordered Multilevel Outcome
%%  ~~function to do ... ~~
}
\description{When several predictors are available, there is often interest in combining a subset of predictors to diagnose disease or predict risk of a clinical outcome, \eqn{D}. In the context of an ordered outcome with \eqn{K} levels, where interest is in predicting \eqn{D = K}, there are multiple ways to select a combination. The traditional approach involves dichotomizing the outcome and using logistic regression to construct the combinations, then selecting a combination based on the estimated AUC for \eqn{D = K} vs. \eqn{D < K} for each fitted combination. An alternative approach, implemented here, constructs the combinations in the same way, but uses both the AUC for \eqn{D = K} vs. \eqn{D < K} and the AUC for \eqn{D = K-1} vs. \eqn{D < K-1}. The combination with the best combined performance is then chosen. This function provides (i) the best combination defined solely by the AUC for \eqn{D=K} vs. \eqn{D < K} and (ii) the best combination defined by both the AUC for \eqn{D = K} vs. \eqn{D < K} and the AUC for \eqn{D = K-1} vs. \eqn{D < K-1}. In the context where \eqn{D} indicates no, mild, or severe disease (\eqn{K}=3), this is equivalent to (i) selecting a combination in terms of its ability to discriminate between individuals with severe vs. no or mild disease and (ii) selecting a combination in terms of its ability to discriminate between individuals with severe vs. no or mild disease and its ability to discriminate between individuals with mild vs. no disease.
}
\usage{
multiselect(data, size=2, Breps=40, nummod=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
      The name of the dataset to be used. An object of class `data.frame' where the first column is the outcome, and the subsequent columns are the predictors. All columns must be numeric. The outcome must be take values 1,...,\eqn{K}, where \eqn{K \geq 3}. Missing observations are not allowed. If the columns of \code{data} are not named, the outcome (first column) will be named "D", and the predictors (subsequent columns) will be named "V1","V2",....
}
\item{size}{
      The size of the combinations. The function considers all possible subsets of the predictors of size \code{size}. Default 2 (all possible pairs).
}
\item{Breps}{
      The number of bootstrap replicates used to estimate the optimism due to resubstitution bias in the AUCs. For each combination, the function estimates the apparent AUCs for each fitted combination. These apparent AUCs are then corrected by substracting the optimism due to resubstitution bias, which is estimated using a bootstrap procedure. Default 40.
}
\item{nummod}{
      The number of predictor combinations to return. Using the optimism-corrected estimate of the AUC for \eqn{D=K} vs. \eqn{D < K}, the function returns the top \code{nummod} predictor combinations. Default 10.
}
}
\details{
For each possible predictor combination of size \code{size}, the function fits the predictor combination using logistic regression comparing outcome \eqn{D = K} to \eqn{D < K}. The apparent AUCs for (a) \eqn{D = K} vs. \eqn{D < K} and (b) \eqn{D = K-1} vs. \eqn{D < K-1} are calculated. A bootstrapping procedure is then used to estimate the optimism due to resubstitution bias in these apparent AUCs. The AUCs are corrected by subtracting the estimated optimism due to resubstitution bias. Two combinations are then selected: the combination with the highest AUC for \eqn{D = K} vs. \eqn{D < K} ("single AUC" approach) and the combination with the best sum of ranks for the AUC for \eqn{D = K} vs. \eqn{D < K} and the AUC for \eqn{D = K-1} vs. \eqn{D < K-1} ("multi-AUC" approach). The selected combinations may be the same for the two approaches. The top \code{nummod} combinations, in terms of the AUC for \eqn{D = K} vs. \eqn{D < K} (corrected for optimism due to resubstitution bias), are also provided.

If more than one combination is "best" in terms of either the AUC for \eqn{D = K} vs. \eqn{D < K} or the sum of ranks for the AUC for \eqn{D = K} vs. \eqn{D < K} and the AUC for \eqn{D = K-1} vs. \eqn{D < K-1} (i.e., in the event of ties) the first combination is returned. The order of the combinations for \eqn{p} candidate predictors is given by \code{combn(1:p, size)}. If ties occur for either (i) the AUC for \eqn{D = K} vs. \eqn{D < K} or (ii) the sum of ranks for the AUC for \eqn{D = K} vs. \eqn{D < K} and the AUC for \eqn{D = K-1} vs. \eqn{D < K-1}, a warning is given.

A given bootstrap sample may not have observations from each of the \eqn{K} outcome levels; if this occurs, a warning is given and the estimated optimism for that bootstrap sample for both the AUC for \eqn{D = K} vs. \eqn{D < K} and the AUC for \eqn{D = K-1} vs. \eqn{D < K-1} will be NA. NAs are removed in the calculation of the mean optimism (used to correct the AUC estimates for resubstitution bias), and the total number of NAs across the \code{Breps} (for either the AUC for \eqn{D = K} vs. \eqn{D < K} or the AUC for \eqn{D = K-1} vs. \eqn{D < K-1}) is indicated by "numNA" in the output.
}
\value{
A list with the following components:
\item{Best.Single }{The best predictor combination as chosen by the "single AUC" approach. The first \code{size} elements give the names of the included predictors (under "Var1","Var2",...), the next is the estimated AUC \eqn{D = K} vs. \eqn{D < K} ("AUC1"), the next is the estimated AUC for \eqn{D = K-1} vs. \eqn{D < K-1} ("AUC2"), the next is number of NAs across the \code{Brep} bootstrap replicates ("numNA"; see 'Details'), and the final \code{size} elements give the estimated coefficients for each of the included predictors ("Coef1","Coef2",...). Both AUCs are corrected for optimism due to resubstitution bias. Recall that if the columns of \code{data} are unnamed, the predictors will be named "V1","V2",....}
\item{Best.Multi }{The best predictor combination as chosen by the "multi-AUC" approach. The elements of \code{Best.Multi} are the same as \code{Best.Single}.}
\item{Ranked.Rslts }{The results for the \code{nummod} best combinations, as ranked by the AUC for \eqn{D = K} vs. \eqn{D < K} (after correcting for optimism due to resubstitution bias). The columns are the same as the elements of \code{Best.Single} and \code{Best.Multi}.
}
}
\references{
{Meisner, A., Parikh, C.R., and Kerr, K.F. (2017). Using multilevel outcomes to construct and select biomarker combinations for single-level prediction. UW Biostatistics Working Paper Series, Working Paper 423.}
}
%% \author{
%%  ~~who you are~~
%% }
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
  \dontshow{
  library(MASS)
  ## shorter example not for users

  set.seed(15)
  p = 4 ## number of predictors
  matX <- matrix(rep(0.3,p*p), nrow=p, ncol=p) ## covariance matrix for the predictors
  diag(matX) <- rep(1,p)

  simD <- apply(rmultinom(400, 1, c(0.6,0.335,0.065)),2,which.max)
  simDord <- simD[order(simD)]
  numobs <- table(simDord)

  simX1 <- mvrnorm(numobs[1], rep(0,p), 2*matX)
  simX2 <- mvrnorm(numobs[2], c(1.5, 1, rep(0.5,(p-2)/2), rep(0.1,(p-2)/2)), 2*matX)
  simX3 <- mvrnorm(numobs[3], c(rep(2,2), rep(0.8,(p-2)/2), rep(0.1,(p-2)/2)), 2*matX)
  simX <- rbind(simX1, simX2, simX3)

  exdata <- data.frame("D"=simDord, simX)

  multiselect(data=exdata, size=2, Breps=2, nummod=3)
  }
  \donttest{
  library(MASS)
  ## example takes ~1 minute to run

  set.seed(15)
  p = 16 ## number of predictors
  matX <- matrix(rep(0.3,p*p), nrow=p, ncol=p) ## covariance matrix for the predictors
  diag(matX) <- rep(1,p)

  simD <- apply(rmultinom(400, 1, c(0.6,0.335,0.065)),2,which.max)
  simDord <- simD[order(simD)]
  numobs <- table(simDord)

  simX1 <- mvrnorm(numobs[1], rep(0,p), 2*matX)
  simX2 <- mvrnorm(numobs[2], c(1.5, 1, rep(0.5,(p-2)/2), rep(0.1,(p-2)/2)), 2*matX)
  simX3 <- mvrnorm(numobs[3], c(rep(2,2), rep(0.8,(p-2)/2), rep(0.1,(p-2)/2)), 2*matX)
  simX <- rbind(simX1, simX2, simX3)

  exdata <- data.frame("D"=simDord, simX)

  multiselect(data=exdata, size=2, Breps=20, nummod=10)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
