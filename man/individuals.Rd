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
'\link{aldmck}', '\link{blackbox}', '\link{blackbox_transpose}'.
}

\examples{
  ### Loads issue scales from the 1980 NES.
  data(Issues1980)
  Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
  Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

  \donttest{ 
  Issues1980_bb <- blackbox(Issues1980, missing=c(0,8,9), verbose=FALSE, 
    dims=3, minscale=8)
  }
  ### 'Issues1980_bb' can be retrieved quickly with: 
  data(Issues1980_bb)

  individuals(Issues1980_bb)
}

\keyword{ multivariate }
