\name{plot.blackbox}
\alias{plot.blackbox}
\title{ Blackbox Coordinate Distribution Plot }
\description{ 
  \code{plot.blackbox} reads an \code{blackbox} object and plots a histogram of the estimated intercepts.
}

\usage{
\method{plot}{blackbox}(x, \ldots)
}

\arguments{
  \item{x}{ an \code{blackbox} output object. }
  \item{...}{ other arguments to \code{hist}.}
}

\value{
    A histogram of the estimated intercepts.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}

}

\seealso{
   '\link{Issues1980}', '\link{summary.blackbox}', '\link{plot.blackbox}'.
}

\examples{

### Loads issue scales from the 1980 NES.
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
result <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
plot(result)

}
\keyword{ multivariate }
