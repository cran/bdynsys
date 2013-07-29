\name{PhasePortraitModelTrajectories}

%\docType{methods}

\alias{phaseportmod}

\title{Phase Potrait with Modeled Trajectories for Specified Entities}

\description{
\code{phaseportmod} creates a phase portrait (velocity plot, vectorfield) based on the two differential equations that were considered to be the overall best models for indicator 1 and indicator 2. The user may specify entities for which he/she wants to obtain the predicted (modeled) trajectories over the phase portrait. The user may extend the period of time for the modeling to get future-predicted trjectories of the entities based on the models. It calls the function \code{preprocess_data}. It requires the packages \code{plm}, \code{pracma}, \code{deSolve} and \code{graphics}.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
\usage{phaseportmod(dataset, yearnr, xv, yv, rangeX, rangeY, param, f, entidx1, 
entidx2, entidx3, entidx4, entidx5, entidx6)}

\arguments{

    \item{dataset}{
        a \code{plm}, a\code{pdata.frame} panel data frame.
        }
    \item{yearnr}{
        number of years for which the user wants to produce the phase portrait
        and obtain the modeled trajectories of the specified entities.
        }
    \item{xv}{
        a reference to variable from the paneldata to be included
        as indicator 1 in the modeling procedure. 
        }
    \item{yv}{
        a reference to variable from the paneldata to be included
        as indicator 2 in the modeling procedure. 
        }
    \item{rangeX}{
        defines the range of the indicator 1. 
        }
    \item{rangeY}{
        defines the range of the indicator 2. 
        }
    \item{param}{
        a vector with parameters from the two models in f. 
        }
    \item{f}{
        a function that contains the models of the two indicators. 
        }
    \item{entidx1}{
        index of an entity in the wide format data of a variable.
        }
    \item{entidx2}{
        index of an entity in the wide format data of a variable.
        }
    \item{entidx3}{
        index of an entity in the wide format data of a variable.
        }
    \item{entidx4}{
        index of an entity in the wide format data of a variable.
        }
    \item{entidx5}{
        index of an entity in the wide format data of a variable.
        }
    \item{entidx6}{
        index of an entity in the wide format data of a variable.
        }  
}

\examples{
## Phase Portrait with Modeles Trajectories over 30 years with the following two models:
## dx  = + 0.0012 /x^2 and dy = + 0.0071 x^3

phaseportmod(datap, 30, datap$logGDP, datap$EmanzV, seq(0, 1, by = 0.1), seq(0, 1, by = 0.1),
param <- c(0.0012, 0.0071), f <- function(t,Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3), 
1, 2, 4, 5, 7, 9)
}

\keyword{graphics}
