\name{FinindexOfModelterms(InternalMethod)}

\alias{findindexM}

\title{Finding Indexes of Polynomial Modelterms}

\description{
\code{findindexM} finds indexs of modelterms to be processed in selectterms if the user wishes to compare Bayes factors of models with two vs. model with three indicators or models with three vs. models with four indicators. The function is internally called in the code \code{bdynsys}{bayesfac}{selectterms}. It returns internally the indexes of modelterms to \code{bdynsys}{bayesfac}{selectterms}.}

\author{Viktoria Spaiser: viktoria.spaiser@iffs.se}
  
%\usage{\method{bdynsys}{bayesfac}{findindexM}(indexM, nterms)}

%\arguments{

%    \item{indexM}{
%        a vector that contains all possible combinations of polynomial terms. 
%        }
%    \item{nterms}{
%        number of all polynomial terms to be possibly used in the modeling. 
%        }
%}

%\examples{
%## Internal call in bdynsys/bayesfac/selectterms
%## paramnumber if the number of model parameters (maximum number of
%## included polynomial terms)
%     wantedarray = sort(c(1, Xterms, Yterms, XYterms)) 
%     Mxy <- c()
%     for (modterm = 1:paramnumber)
%     {
%       comb <- combs(wantedarray, modterm)
%       indexset = c()
%       for (t in 1:nrow(comb))
%       {
%         indexM <- comb[t, 1:modterm]
%         m <- findindexM(indexM, nterms)
%         indexset <- c(indexset, m)
%       } 
%       ind = which.min(SqErX[modterm, indexset])
%       Mxy[modterm] <- indexset[ind]      
%     } 