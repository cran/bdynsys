\name{BayesianModelFitting(InternalMethod)}

\alias{polyfitbayes}

\title{Bayesian Model Fitting with Polynomial Terms als Predictors and Changes as Dependent Variable}

\description{
\code{polyfitbayes} is a modeling function internally called in the code \code{bdynsys}{bayesfac}. It computes the values for the independed variables, based on polynomial terms and a coefficient parameter space. The function then performs a Monte Carlo simulation to compute/integrate the log Bayes Factor. It returns internally the log Bayes factors for the best models and the modelindexes of those best models to \code{bdynsys}{bayesfac}.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
%\usage{\method{bdynsys}{bayesfac}{polyfitbayes}(indicatornumber, parameternumber, x, y, ch, selected_terms, z, v)}
%
%\arguments{
%
%    \item{indicatornumber}{
%        an integer number indicating number of indicators,
%        to be included in the modeling procedure 
%        }
%    \item{parameternumber}{
%        an integer number indicating number of modelparameters,
%        this is the maximum number of polynomial terms included. 
%        }
%   \item{x}{
%        a vector data object passed on from \code{bayesfac},  
%        which is referring to indicator 1 in the modeling procedure. 
%        }
%   \item{y}{
%        a vector data object passed on from \code{bayesfac},  
%        which is referring to indicator 2 in the modeling procedure. 
%        }
%   \item{ch}{
%        a vector data object passed on from \code{bayesfac},
%        which is referring to changes in a respective indicator at each timestep. 
%        }
%   \item{selected_terms}{
%        a vector that contains the respective combination of polynimal terms, passed on from 
%        \code{bayesfac}. 
%        }
%   \item{z}{
%        a vector data object passed on from \code{bayesfac},  
%        which is referring to indicator 3 in the modeling procedure. 
%        }
%   \item{v}{
%        a vector data object passed on from \code{bayesfac},  
%        which is referring to indicator 4 in the modeling procedure. 
%        }   
%}
%
%\examples{
%## Internal call in bdynsys/bayesfac with two indicators (indnr) 
%## here for dx/dt models
%## paramnumber if the number of model parameters (maximum number of
%## included polynomial terms)
%BayesFacX <- c()
%parambestx <- c()   
%Mx <- Mx_allterms
%count = 1   
%for (j in 1:paramnumber)
%{
%  tmp = polyfitbayes(indnr, x, y, chx, M1[count:count+j-1]) 
%  bestmbf <- tmp[[1]]
%  indexbestm <- tmp[[2]]
%      
%  BayesFacX[j] <- bestm     
%  parambestx[j] <- indexbestm
%  count = count + j;
%}
%
%\keyword{stats}
