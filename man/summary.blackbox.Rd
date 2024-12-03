\name{summary.blackbox}
\alias{summary.blackbox}
\title{ Blackbox Summary }
\description{ 
  \code{summary.blackbox} reads an \code{blackbox} object and prints a summary.
}

\usage{
  \method{summary}{blackbox}(object, \ldots)             
}

\arguments{
  \item{object}{ a \code{blackbox} output object. }
  \item{...}{ further arguments to \code{print}. }
}

\value{
    A summary of a \code{blackbox} object. For each dimension, reports all
    stimuli with coordinates, individuals used for scaling, and fit. Also
    summarizes number of rows, columns, total data entries, number of missing
    entries, percent missing data, and sum of squares.
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

  \donttest{ 
  Issues1980_bb <- blackbox(Issues1980, missing=c(0,8,9), verbose=FALSE, 
    dims=3, minscale=8)
  }
  ### 'Issues1980_bb' can be retrieved quickly with: 
  data(Issues1980_bb)

  summary(Issues1980_bb)
}

\keyword{ multivariate }
