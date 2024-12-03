\name{predict.blackbox}
\alias{predict.blackbox}
\title{ Predict method of blackbox objects }
\description{ 
  \code{predict.blackbox} reads an \code{blackbox} object and uses the estimates to generate a matrix of predicted values.
}

\usage{
  \method{predict}{blackbox}(object, dims=1, ...)             
}

\arguments{
  \item{object}{ A \code{blackbox} output object. }
  \item{dims}{ Number of dimensions used in prediction. Must be equal to or less than number of dimensions used in estimation. }
  \item{...}{ Ignored. }
}

\value{
    A matrix of predicted values generated from the parameters estimated from a \code{blackbox} object.
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
David A. Armstrong II, Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2021. \emph{Analyzing Spatial Models of Choice and Judgment}. 2nd ed. Statistics in the Social and Behavioral Sciences Series. Boca Raton, FL: Chapman & Hall/CRC. doi: 10.1201/9781315197609

Keith T. Poole, Jeffrey B. Lewis, Howard Rosenthal, James Lo, and Royce Carroll. 2016. ``Recovering a Basic Space from Issue Scales in R.'' \emph{Journal of Statistical Software} 69(7): 1-21. doi:10.18637/jss.v069.i07

Keith T. Poole. 1998. ``Recovering a Basic Space From a Set of Issue Scales.'' \emph{American Journal of Political Science} 42(3): 954-993. doi: 10.2307/2991737
}


\seealso{
'\link{blackbox}', '\link{Issues1980}'
}

\examples{
  ### Loads issue scales from the 1980 ANES.
  data(Issues1980)
  Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
  Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

  ### Estimate blackbox object from example and call predict function
  \donttest{ 
  Issues1980_bb <- blackbox(Issues1980, missing=c(0,8,9), verbose=FALSE, 
    dims=3, minscale=8)
  }
  ### 'Issues1980_bb' can be retrieved quickly with: 
  data(Issues1980_bb)
  prediction <- predict.blackbox(Issues1980_bb, dims=3)

  ### Examine predicted vs. observed values for first 10 respondents
  ### Note that 4th and 6th respondents are NA because of missing data
  Issues1980[1:10,]
  prediction[1:10,]

  ### Check correlation across all predicted vs. observed, excluding missing values
  prediction[which(Issues1980 \%in\% c(0,8,9))] <- NA
  cor(as.numeric(prediction), as.numeric(Issues1980), use="pairwise.complete")
}

\keyword{ multivariate }
