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

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{aldmck}', '\link{LC1980}'
}

\examples{
## Estimate an aldmck object from example and call predict function
data(LC1980)
result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
prediction <- predict.aldmck(result)

## Examine predicted vs. observed values for first 10 respondents
## Note some observations are NA'd in prediction matrix from caliper
## First column of LC1980 are self-placements, which are excluded
LC1980[1:10,-1]
prediction[1:10,]

## Check correlation across all predicted vs. observed, excluding missing values
prediction[which(LC1980[,-1] \%in\% c(0,8,9))] <- NA
cor(as.numeric(prediction), as.numeric(LC1980[,-1]), use="pairwise.complete")
}
\keyword{ multivariate }
