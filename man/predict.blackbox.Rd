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

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{blackbox}', '\link{Issues1980}'
}

\examples{
## Estimate blackbox object from example and call predict function
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

### This command conducts estimates, which we instead load using data()
# Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
data(Issues1980_bb)
prediction <- predict.blackbox(Issues1980_bb,dims=3)

## Examine predicted vs. observed values for first 10 respondents
## Note that 4th and 6th respondents are NA because of missing data
Issues1980[1:10,]
prediction[1:10,]

## Check correlation across all predicted vs. observed, excluding missing values
prediction[which(Issues1980 \%in\% c(0,8,9))] <- NA
cor(as.numeric(prediction), as.numeric(Issues1980), use="pairwise.complete")
}
\keyword{ multivariate }
