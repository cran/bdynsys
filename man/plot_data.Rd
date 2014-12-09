\name{plot_data}

%\docType{methods}

\alias{plot_data}

\title{Plotting Data in Phase Plane with Data Trajectories for Specified Entities}

\description{
\code{plot_data} plots data form indicator 1 and indicator 2 in phase plane. The user may specify entities for which he/she wants to plot the data trajectories over the plotted data. It calls the function \code{preprocess_data}. It requires the packages \code{plm} and \code{graphics}.}

\author{Viktoria Spaiser: viktoria.sp@web.de}
  
\usage{plot_data(dataset, xvar, yvar, rangeX, rangeY, entidx1, entidx2, entidx3, 
entidx4, entidx5, entidx6)}

\arguments{

    \item{dataset}{
        a \code{plm} \code{pdata.frame} panel data frame.
        }
    \item{xvar}{
        a reference to variable from the paneldata to be included
        as indicator 1 in the modeling procedure. 
        }
    \item{yvar}{
        a reference to variable from the paneldata to be included
        as indicator 2 in the modeling procedure. 
        }
    \item{rangeX}{
        defines the range of the indicator 1. 
        }
    \item{rangeY}{
        defines the range of the indicator 2. 
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
## Plot Data in a Phase Plane with Data Trajectories for Specified Entities

plot_data(datap, datap$logGDP, datap$EmanzV, seq(0, 12, by = 0.5), seq(0, 1, by = 0.1), 
1, 2, 3, 4, 5, 6)
}

\keyword{graphics}