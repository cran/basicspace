\name{predict.aldmck}
\alias{predict.aldmck}
\title{ Predict method of aldmck objects }
\description{ 
  \code{predict.aldmck} reads an \code{aldmck} object and uses the estimates to generate a matrix of predicted values.
}

\usage{
  \method{predict}{aldmck}(object, caliper=0.2, ...)             
}

\arguments{
  \item{object}{ A \code{aldmck} output object. }
  \item{caliper}{ Caliper tolerance. Any individuals with estimated weights lower than this value are NA'd out for prediction. Since predictions are made by dividing observed values by estimating weights, very small weights will grossly inflate the magnitude of predicted values and lead to extreme predictions.}
  \item{...}{ Ignored.}
}

\value{
    A matrix of predicted values generated from the parameters estimated from a \code{aldmck} object.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}

Christopher Hare \email{cdhare@ucdavis.edu}
}

\references{
John H. Aldrich and Richard D. McKelvey. 1977. ``A Method of Scaling with Applications to the 1968 and 1972 Presidential Elections.'' \emph{American Political Science Review} 71(1): 111-130. doi: 10.2307/1956957

David A. Armstrong II, Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2021. \emph{Analyzing Spatial Models of Choice and Judgment}. 2nd ed. Statistics in the Social and Behavioral Sciences Series. Boca Raton, FL: Chapman & Hall/CRC. doi: 10.1201/9781315197609

Thomas R. Palfrey and Keith T. Poole. 1987. ``The Relationship between Information, Ideology, and Voting Behavior.'' \emph{American Journal of Political Science} 31(3): 511-530. doi: 10.2307/2111281

Keith T. Poole, Jeffrey B. Lewis, Howard Rosenthal, James Lo, and Royce Carroll. 2016. ``Recovering a Basic Space from Issue Scales in R.'' \emph{Journal of Statistical Software} 69(7): 1-21. doi:10.18637/jss.v069.i07

Keith T. Poole. 1998. ``Recovering a Basic Space From a Set of Issue Scales.'' \emph{American Journal of Political Science} 42(3): 954-993. doi: 10.2307/2991737
}

\seealso{
'\link{aldmck}', '\link{LC1980}'
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)

  ### Estimate an aldmck object from example and call predict function
  result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), verbose=TRUE)
  prediction <- predict.aldmck(result)

  ### Examine predicted vs. observed values for first 10 respondents
  ### Note some observations are NA'd in prediction matrix from caliper
  ### First column of LC1980 are self-placements, which are excluded
  LC1980[1:10,-1]
  prediction[1:10,]

  ### Check correlation across all predicted vs. observed, excluding missing values
  prediction[which(LC1980[,-1] \%in\% c(0,8,9))] <- NA
  cor(as.numeric(prediction), as.numeric(LC1980[,-1]), use="pairwise.complete")
}

\keyword{ multivariate }
