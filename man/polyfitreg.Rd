\name{PolynomialTermsRegression(InternalMethod)}

\alias{polyfitreg}
\alias{regression}

\title{Least Square Regression with Polynomial Terms als Predictors and Changes as Dependent Variable}

\description{
\code{polyfitreg} is a modeling function internally called in the code \code{bdynsys}{dysymod} as well as in \code{bdynsys}{dysymod}{bestfitmod}. It computes the values for the independed variables, based on polynomial terms and computes all possible regression models with all possible combinations of these independ variables (all possible combinations of polynimal terms) and the changes in the dependent variables which are to be predicted. The function estimates the Betas, computes the Sum Square of Error of all possible models as well as their Log Likelihood and R-squares. It returns internally the Sum Square of Error matrices to \code{bdynsys}{dysymod}{bestfit} and Betas to \code{bdynsys}{dysymod}.Log Likelihood and R-square are printed.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
%\usage{\method{bdynsys}{dysymod}{polyfitreg}(indicatornumber, parameternumber, x, y, ch, selected_terms, z, v)}
%  
%\arguments{
%
%    \item{indicatornumber}{
%        an integer number indicating number of indicators,
%        to be included in the modeling procedure 
%       }
%    \item{parameternumber}{
%        an integer number indicating number of modelparameters,
%        this is the maximum number of polynomial terms included. 
%        }
%   \item{x}{
%        a vector data object passed on from \code{dysymod} and from \code{dysymod}{bestfit},  
%        which is referring to indicator 1 in the modeling procedure. 
%        }
%   \item{y}{
%        a vector data object passed on from \code{dysymod} and from \code{dysymod}{bestfit},  
%        which is referring to indicator 2 in the modeling procedure. 
%        }
%   \item{ch}{
%        a vector data object passed on from \code{dysymod} and from \code{dysymod}{bestfit},
%        which is referring to changes in a respective indicator at each timestep. 
%        }
%   \item{selected_terms}{
%        a vector that contains the respective combination of polynimal terms, passed on from 
%        \code{dysymod} and from \code{dysymod}{bestfit}. 
%        }
%   \item{z}{
%        a vector data object passed on from \code{dysymod} and from \code{dysymod}{bestfit},  
%        which is referring to indicator 3 in the modeling procedure. 
%        }
%   \item{v}{
%        a vector data object passed on from \code{dysymod} and from \code{dysymod}{bestfit},  
%        which is referring to indicator 4 in the modeling procedure. 
%        }   
%}
%
%\examples{
%## Internal call in bdynsys/dysymod/bestfitmod with two indicators (indnr) 
%## paramnumber if the number of model parameters (maximum number of
%## included polynomial terms)
%## nterms is the number of all possible predifined polynimal terms
%for (i in 1:parameternumber)
%    {
%      M <- combs(1:nterms, i)     
%      for (j in 1:nrow(M)) 
%      {
%        tmp = polyfitreg(indnr, x, y, ch, M[j,])
%        Betas <- tmp[[1]]
%        SSqEr <- tmp[[2]]        
%      }  
%    }
%}
%
%\keyword{stats
%models}
