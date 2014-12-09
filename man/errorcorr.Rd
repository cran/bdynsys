\name{errorcorr}

%\docType{methods}

\alias{errorcorr}

\title{Controlling Error Correlations in Models with Panel Data}

\description{
\code{errorcorr} is an additional function in the bdynsys package. It calls functions \code{preprocess_data}. It computes the errors in the models (differential equations) and the covariances of the errors. It then uses the covarinaces to re-estimate the Betas of the models. The covarinace matrix is printed and the re-estimated Betas saved in a file. It requires the package \code{MASS}}

\author{
Viktoria Spaiser: viktoria.sp@web.de}

\usage{errorcorr(dataset, indnr, x, y, f, xterms, yterms, nrterms, z, zterms, v, vterms)}
	
\arguments{

    \item{dataset}{
        a \code{plm} \code{pdata.frame} panel data frame.
        }
    \item{indnr}{
        an integer number indicating number of indicators,
        to be included in the modeling procedure 
        }
   \item{x}{
        a reference to variable from the paneldata to be included
        as indicator 1 in the modeling procedure. 
        }
   \item{y}{
        a reference to variable from the paneldata to be included
        as indicator 2 in the modeling procedure. 
        }
   \item{f}{
        a function that contains the models of the indicators. 
        }
   \item{xterms}{
        a vector that contains the terms from the model dx/dt. 
        }
   \item{yterms}{
        a vector that contains the terms from the model dy/dt. 
        }
   \item{nrterms}{
        total number of in all equations, e.g. sum of terms in
        equation for dx/dt and terms in equation for dy/dt, if 
        the number of variables is two. 
        }
   \item{z}{
        a reference to variable from the paneldata to be included
        as indicator 3 in the modeling procedure. 
        }
   \item{zterms}{
        a vector that contains the terms from the model dz/dt. 
        }
   \item{v}{
        a reference to variable from the paneldata to be included
        as indicator 4 in the modeling procedure. 
        } 
   \item{vterms}{
        a vector that contains the terms from the model dv/dt. 
        }
}

\examples{
## Controlling Error Correlations with two indicators and with the following two models:
## dx/dt  = + 0.0012 /x^2 and dy/dt = + 0.0071 x^3

errorcorr(datap, 2, datap$logGDP, datap$EmanzV, 
f <- function(Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3), c(11), c(14), 2)
}

\keyword{stats}
