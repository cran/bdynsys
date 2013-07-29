\name{BayesianFactorComputation(InternalMethod)}

\alias{bayesfac}
\alias{modelselection}

\title{Bayesian Model Selection with Bayesian Factor Computation}

\description{
\code{bayesfac} is the core code of the bayesian model selection in the bdynsys package. It is internally called in the main code \code{bdynsys}. It calls the functions \code{polyfitbayes}, \code{selectterms} and \code{findindexM}. It uses the Sum Square Errors from the previous model fitting step in \code{dysymod} to identify the best models for each number of modelparameters and computes then the Bayesian factor using a Monte Carlo Simulation for integration for these models, which may be then compared to find the overall best model. It prints the Bayesian factors for the best models for each number of modelparameters.It requires the package \code{caTools}.}

\author{Viktoria Spaiser, viktoria.spaiser@iffs.se}