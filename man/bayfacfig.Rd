\name{BayesFactorsPlot}

\alias{bayfacfig}

\title{Bayes Factors Plot}

\description{
\code{bayfacfig} plot the Bayes Factors to visualy display the selected and the overall best model. It requires the packages \code{Hmisc} and \code{graphics}.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
\usage{bayfacfig(indnr, modelnr, BF, markmod)}

\arguments{

    \item{indnr}{
        an integer number indicating number of indicators,
        to be included in the modeling procedure 
        }
    \item{modelnr}{
        number of models for which Bayes factors were obtained,
        and are to be compared now in the plot
        }
    \item{BF}{
        double numbers giving the Bayes factors for the models,
        to be compared, stored in a vector.  
        }
    \item{markmod}{
        the number of the model (index of the Bayes factor in the vector),
        which is the overall best. 
        }  
}

\examples{
## Bayes Factor plot for five models with two indicators to be compared in the 
## bayesian fitting, with forth model being highlighted as the overall best model 
bayfacfig(2, 5, c(-5.4534, -5.3955, -5.235, -4.99948, -5.321), 4)
}

\keyword{graphics}
