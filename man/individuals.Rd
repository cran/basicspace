\name{individuals}
\alias{individuals}
\title{ Extraction function for scaled individuals}
\description{ 
  \code{individuals} is a convenience function to extract the individual/respondent parameters from an \code{aldmck}, \code{blackbox}, or \code{blackbt} object.
}

\usage{
  individuals(object)
}

\arguments{
  \item{object}{ an \code{aldmck}, \code{blackbox}, or \code{blackbt} output object. }
}

\value{
    The individual parameters of the estimated output, which can also be recovered as
    \code{object$individuals} (for \code{blackbox} or \code{blackbt} objects) or
    \code{object$respondents} (for \code{aldmck} objects).  Please refer to the
    documentation of \code{aldmck}, \code{blackbox}, or \code{blackbox_transpose}
    for specifics.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{aldmck}', '\link{blackbox}', '\link{blackbox_transpose}'.
}

\examples{
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

### This command conducts estimates, which we instead load using data()
# Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
data(Issues1980_bb)

individuals(Issues1980_bb)
}
\keyword{ multivariate }
