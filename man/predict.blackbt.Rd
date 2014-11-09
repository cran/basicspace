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

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{blackbox_transpose}', '\link{LC1980}', '\link{LC1980_bbt}'
}

\examples{
## Estimate blackbt object from example and call predict function
data(LC1980)
data(LC1980_bbt)
prediction <- predict.blackbt(LC1980_bbt, dims=2)

## Examine predicted vs. observed values for first 10 respondents
## First column of LC1980 are self-placements, which are excluded
LC1980[1:10,-1]
prediction[1:10,]

## Check correlation across all predicted vs. observed, excluding missing values
prediction[which(LC1980[,-1] \%in\% c(0,8,9))] <- NA
cor(as.numeric(prediction), as.numeric(LC1980[,-1]), use="pairwise.complete")
}
\keyword{ multivariate }
