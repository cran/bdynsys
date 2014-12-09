\name{dysymod}

\alias{dysymod}

\title{Modeling Dynamic Complex Systems with Ordinary Differential Equations (Internal Function)}

\description{
\code{dysymod} is the core code of the dynamic system modeling in the bdynsys package. It is internally called in the main code \code{bdynsys}. It calls the functions \code{polyfitreg} and \code{bestfitmod}. It defines the polynomial terms, that are used to build the models and creates all possible model combinations and fits them. The best three models for each number of modelparameters are then selected. It prints the three selected models for each number of modelparameters and their Log Likelihood and R-squares. It requires the package \code{caTools}.}

\author{Viktoria Spaiser: viktoria.sp@web.de}
 
%\usage{\method{bdynsys}{dysymod}(indicatornumber, parameternumber, x, y, chx, chy, mx, my, z, chz, mz, v, chv, mv)}
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
%        a return object from \code{preprocess_data}, which is referring
%        to indicator 1 in the modeling procedure. 
%        }
%   \item{y}{
%        a return object from \code{preprocess_data}, which is referring
%        to indicator 2 in the modeling procedure. 
%        }
%   \item{chx}{
%        a return object from \code{preprocess_data}, which is referring
%        to changes in indicator 1 at each timestep. 
%        }
%   \item{chy}{
%        a return object from \code{preprocess_data}, which is referring
%        to changes in indicator 2 at each timestep. 
%        }
%   \item{mx}{
%        a return object from \code{preprocess_data}, which is referring
%        to the mean of indicator 1. 
%        }
%   \item{my}{
%        a return object from \code{preprocess_data}, which is referring
%        to the mean of indicator 2. 
%        }
%   \item{z}{
%        a return object from \code{preprocess_data}, which is referring
%        to indicator 3 in the modeling procedure. 
%        }
%   \item{chz}{
%        a return object from \code{preprocess_data}, which is referring
%        changes in indicator 3 at each timestep. 
%        }
%   \item{mz}{
%        a return object from \code{preprocess_data}, which is referring
%        to the mean of indicator 3. 
%        }
%  \item{v}{
%        a return object from \code{preprocess_data}, which is referring
%        to indicator 4 in the modeling procedure. 
%        } 
%  \item{chv}{
%        a return object from \code{preprocess_data}, which is referring
%        changes in indicator 4 at each timestep. 
%        } 
%  \item{mv}{
%        a return object from \code{preprocess_data}, which is referring
%        to the mean of indicator 4. 
%        }
%}
%
%\examples{
%## Internal Call in bdynsys with two indicators  
%results <- dysymod(indnr, paramnr, procdata$xs, procdata$ys, procdata$chXs, procdata$chYs, 
%                    procdata$mx, procdata$my)
%
%## Internal Call in bdynsys with three indicators 
%results <- dysymod(indnr, paramnr, procdata$xs, procdata$ys, procdata$chXs, procdata$chYs, 
%                    procdata$mx, procdata$my, procdata$zs, procdata$chZs, procdata$mz)
%
%## Internal Call in bdynsys with four indicators
%results <- dysymod(indnr, paramnr, procdata$xs, procdata$ys, procdata$chXs, procdata$chYs, 
%                    procdata$mx, procdata$my, procdata$zs, procdata$chZs, procdata$mz,
%                    procdata$vs, procdata$chVs, procdata$mv)
%}
%
%\keyword{stats
%models
%nonlinear
%ts
%multivariate}
