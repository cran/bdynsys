\name{SelectTerms(InternalMethod)}

\alias{selectterms}

\title{Deriving Best Models and their Polynomial Terms for Bayesian Model Fitting}

\description{
\code{selectterms} selects/indexes terms which are part of the models to be fitted in the Bayesian fitting procedure. The function is internally called in the code \code{bdynsys}{bayesfac}. The model selections is based on the Sum of Square Errors which were computed in the previous model fitting step in \code{bdynsys}{dysymod}. It returns internally the indexes of the best models polynomial terms to \code{bdynsys}{bayesfac}.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
%\usage{\method{bdynsys}{bayesfac}{selectterms}(indicatornumber, parameternumber, SSEx, SSEy, SSEz, SSEv)}
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
%    \item{SSEx}{
%        an object (matrix) passed on from \code{bayesfac} giving the Sum of Square
%        Errors for all possible models for dx/dt.  
%        }
%    \item{SSEy}{
%        an object (matrix) passed on from \code{bayesfac} giving the Sum of Square
%        Errors for all possible models for dy/dt.  
%        }
%   \item{SSEz}{
%        an object (matrix) passed on from \code{bayesfac} giving the Sum of Square
%        Errors for all possible models for dz/dt.  
%        }
%   \item{SSEv}{
%        an object (matrix) passed on from \code{bayesfac} giving the Sum of Square
%        Errors for all possible models for dv/dt.  
%        }   
%}
%
%\examples{
%## Internal call in bdynsys/bayesfac with two indicators (indnr) 
%tmp = selectterms(indnr, paramnr, SSEx, SSEy) 
%Mx <- tmp[[1]]
%My <- tmp[[2]]
%
%\keyword{stats}