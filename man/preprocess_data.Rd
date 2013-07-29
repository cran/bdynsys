\name{DataProcessing(InternalMethod)}

\alias{preprocess_data}

\title{Data Processing for Bayesian Dynamic System Modeling}

\description{
\code{preprocess_data} is a core code in the bdynsys package to prepare the data for the bayesian dynamic system modeling. It is internally called in the main code \code{bdynsys}. It creates wide format matrixes for the indicators based on time and entity-id. It also reshapes the indicators into single columns ordered by time, computes the changes in the indicators at each timestep, reshapes them also into single colums ordered by time and stores them in variables. It furthermore removes missing values and scales the variables with their means. It returnes the processed variables to \code{bdynsys} and that way makes them availabe for usage in \code{bdynsys}{dysymod} and \code{bdynsys}{bayesfac}.It requires the package \code{plm} and \code{matrixStats}.}

\author{
Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
%\usage{\method{bdynsys}{preprocess_data}(indicatornumber, x, y, z, v)}
%  
%\arguments{
%
%    \item{indicatornumber}{
%        an integer number indicating number of indicators,
%        to be included in the modeling procedure 
%        }
%   \item{x}{
%        a data frame object from \code{bdynsys}, which is referring
%        to indicator 1 in the modeling procedure. 
%        }
%   \item{y}{
%        a data frame object from \code{bdynsys}, which is referring
%        to indicator 2 in the modeling procedure. 
%        }
%   \item{z}{
%        a data frame object from \code{bdynsys}, which is referring
%        to indicator 3 in the modeling procedure. 
%        }
%   \item{v}{
%        a data frame object from \code{bdynsys}, which is referring
%        to indicator 4 in the modeling procedure. 
%        }   
%}
%
%\examples{
%## Internal Call in bdynsys with two indicators  
%procdata <- preprocess_data(indnr, x, y)
%
%## Internal Call in bdynsys with three indicators 
%procdata <- preprocess_data(indnr, x, y, z)
%
%## Internal Call in bdynsys with four indicators
%procdata <- preprocess_data(indnr, x, y, z, v)
%}
%
%\keyword{data}
