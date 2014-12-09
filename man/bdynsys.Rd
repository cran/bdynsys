\name{bdynsys}

\docType{package}

\alias{bdynsys}

\title{Bayesian Dynamic System Modeling}

\description{
\code{bdynsys} is a package for panel/longitudinal data that combines methods to model 
changes in up to four indicators over time as a function of the indicators itself and up to four predictors using ordinary differential equations (ODEs) with polynomial terms that allow to model complex and nonlinear effects. A Bayesian model selection approach is implemented. The package provides also visualisation tools to plot phase portraits of the dynamic system, showing the complex co-evolution of two indicators over time with the possibility to highlight trajectories for specified entities (e.g. countries, individuals). \code{bdynsys} is also the
name of the main function in the \code{bdynsys} package, that performs the bayesian dynamic
systems modeling.}

\references{
       Ranganathan, S./Spaiser, V./Mann, R.P./Sumpter, D.J.T. (2014) 
       \emph{Bayesian Dynamical Systems Modelling in the Social Sciences}.
       PLoS ONE, 9(1):e86468.
     }

\author{Viktoria Spaiser: viktoria.sp@web.de}

\usage{bdynsys(dataset, indnr, paramnr, x, y, z, v)
}
  
\arguments{

    \item{dataset}{
        a \code{plm} \code{pdata.frame} panel data frame.
        }
    \item{indnr}{
        an integer number indicating number of indicators,
        to be included in the modeling procedure 
        }
    \item{paramnr}{
        an integer number indicating number of modelparameters,
        this is the maximum number of polynomial terms included. 
        }
   \item{x}{
        a reference to variable from the paneldata to be included
        as indicator 1 in the modeling procedure. 
        }
   \item{y}{
        a reference to variable from the paneldata to be included
        as indicator 2 in the modeling procedure. 
        }
   \item{z}{
        a reference to variable from the paneldata to be included
        as indicator 3 in the modeling procedure. 
        }
   \item{v}{
        a reference to variable from the paneldata to be included
        as indicator 4 in the modeling procedure. 
        } 
}

\examples{
## Bayesian Dynamic System Modeling with two variables and one modelparameter

bdynsys(datap, 2, 1, datap$logGDP, datap$EmanzV)
}

\keyword{stats}
