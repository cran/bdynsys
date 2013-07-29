\name{FittingDynamicSystemModels(InternalMethod)}

\alias{bestfitmod}
\alias{modelfitting}

\title{Fitting of Dynamic System Models}

\description{
\code{bestfitmod} is a model fitting function internally called in the code \code{bdynsys}{dysymod}. It calls the functions \code{polyfitreg}. It creates the testing and training datasubsets for model fitting, passes them on to \code{polyfitreg} and receives back the Sum Square of Errors which are then stored in matrix for each indicator. It returns internally the Sum Square of Error matrices to \code{bdynsys}{dysymod}. It requires the package \code{caTools}.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
%\usage{\method{bdynsys}{dysymod}{bestfitmod}(indicatornumber, parameternumber, x, y, ch, z, v)}
  
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
%        a vector data object passed on from \code{dysymod}, which is referring
%        to indicator 1 in the modeling procedure. 
%        }
%   \item{y}{
%        a vector data object passed on from \code{dysymod}, which is referring
%        to indicator 2 in the modeling procedure. 
%        }
%   \item{ch}{
%        a vector data object passed on from \code{dysymod}, which is referring
%        changes in a respective indicator at each timestep. 
%        }
%   \item{z}{
%        a vector data object passed on from \code{dysymod}, which is referring
%        to indicator 3 in the modeling procedure. 
%        }
%  \item{v}{
%        a vector data object passed on from \code{dysymod}, which is referring
%        to indicator 4 in the modeling procedure. 
%        }   
%}

%\examples{
%## Internal Call in bdynsys/dysymod for the two indicators   
%    SEtestx <- (bestfitmod(indnr, paramnr, var1, var2, chVar1))
%    SEtesty <- (bestfitmod(indnr, paramnr, var1, var2, chVar2))
%}

%\keyword{stats
%models}
