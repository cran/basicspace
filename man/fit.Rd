\name{fit}
\alias{fit}
\title{ Extraction function for fit of scaling model}
\description{ 
  \code{fit} is a convenience function to extract the model fit statistics from an \code{aldmck}, \code{blackbox}, or \code{blackbt} object.
}

\usage{
  fit(object)
}

\arguments{
  \item{object}{ an \code{aldmck}, \code{blackbox}, or \code{blackbt} output object. }
}

\value{
    The model fit statistics of the estimated output, which can also be recovered as
    \code{object$fits} (for \code{blackbox} or \code{blackbt} objects) or
    \code{object$AMfit} (for \code{aldmck} objects).  Please refer to the
    documentation of \code{aldmck}, \code{blackbox}, or \code{blackbox_transpose}
    for specifics.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{aldmck}', '\link{blackbox}', '\link{blackbox_transpose}'.
}

\examples{
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
result <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
fit(result)
}
\keyword{ multivariate }
