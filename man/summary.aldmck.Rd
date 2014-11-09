\name{summary.aldmck}
\alias{summary.aldmck}
\title{ Aldrich-McKelvey Summary }
\description{ 
  \code{summary.aldmck} reads an \code{aldmck} object and prints a summary.
}

\usage{
\method{summary}{aldmck}(object, \ldots)              
}

\arguments{
  \item{object}{ an \code{aldmck} output object. }
  \item{...}{ further arguments to \code{print}. }
}

\value{
    A summary of an \code{aldmck} object. Reports number of stimuli, respondents
    scaled, number of respondents with positive and negative weights, R-squared,
    Reudction of normalized variance of perceptions, and stimuli locations.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{aldmck}', '\link{LC1980}', '\link{plot.aldmck}', '\link{plot.cdf}'.
}

\examples{
data(LC1980)
result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
summary(result)
par(ask=TRUE)
plot.AM(result,xlim=c(-1.5,1.5))
plot.cdf(result)
}
\keyword{ multivariate }
