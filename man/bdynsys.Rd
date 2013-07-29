\name{BayesianDynamicSystemModeling}

\docType{package}

\alias{bdynsys}
\alias{panel}
\alias{longitudinal}
\alias{bayesian}
\alias{dynamicSystems}
\alias{complex}
\alias{nonlinear}

\title{Bayesian Dynamic System Modeling}

\description{
\code{bdynsys} is a package for panel/longitudinal data that combines methods to model 
changes in up to four indicators over times as a function of the indicators itself and up to three predictors using ordinary differential equations (ODEs) with polynomial terms that allow to model complex and nonlinear effects. A Bayesian model selection approach is implemented. The package provides also visualisation tools to plot phase portraits of the dynamic system, showing the complex co-evolution of two indicators over time with the possibility to highlight trajectories for specified entities (e.g. countries, individuals). \code{bdynsys} is also the
name of the main function in the \code{bdynsys} package, that performs the bayesian dynamic
systems modeling.}

\references{
       Ranganathan, S./Spaiser, V./Sumpter, D.J.T. (2013) 
       \emph{A Bayesian approach to modeling dynamical systems in the social sciences}.
       SIMULTECH 2013 Proceedings, Reykjavik.
     }

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}

\usage{bdynsys(dataset, indnr, paramnr, x, y, z, v)
}
  
\arguments{

    \item{dataset}{
        a \code{plm}, a\code{pdata.frame} panel data frame.
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
## Bayesian Dynamic System Modeling with two variables and three modelparameters
## dx  = + 0.0012 /x^2 and dy = + 0.0071 x^3

bdynsys(datap, 2, 3, datap$logGDP, datap$EmanzV)
}

\keyword{stats}
