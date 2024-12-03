\name{predict.blackbt}
\alias{predict.blackbt}
\title{ Predict method of blackbt objects }
\description{ 
  \code{predict.blackbt} reads an \code{blackbt} object and uses the estimates to generate a matrix of predicted values.
}

\usage{
  \method{predict}{blackbt}(object, dims=1, ...)             
}

\arguments{
  \item{object}{ A \code{blackbox} output object. }
  \item{dims}{ Number of dimensions used in prediction. Must be equal to or less than number of dimensions used in estimation. }
  \item{...}{ Ignored. }
}

\value{
    A matrix of predicted values generated from the parameters estimated from a \code{blackbt} object.
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
'\link{blackbox_transpose}', '\link{LC1980}', '\link{LC1980_bbt}'
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  LCdat <- LC1980[,-1]	#Dump the column of self-placements

  ### Estimate blackbt object from example and call predict function
  \donttest{ 
  LC1980_bbt <- blackbox_transpose(LCdat, missing=c(0,8,9), dims=3, 
    minscale=5, verbose=TRUE)
  }
  ### 'LC1980_bbt' can be retrieved quickly with: 
  data(LC1980_bbt)
  prediction <- predict.blackbt(LC1980_bbt, dims=2)

  ### Examine predicted vs. observed values for first 10 respondents
  ### First column of LC1980 are self-placements, which are excluded
  LC1980[1:10,-1]
  prediction[1:10,]

  ### Check correlation across all predicted vs. observed, excluding missing values
  prediction[which(LC1980[,-1] \%in\% c(0,8,9))] <- NA
  cor(as.numeric(prediction), as.numeric(LC1980[,-1]), use="pairwise.complete")
}

\keyword{ multivariate }
