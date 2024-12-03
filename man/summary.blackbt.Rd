\name{summary.blackbt}
\alias{summary.blackbt}
\title{ Blackbox-Transpose Summary }
\description{ 
  \code{summary.blackbt} reads an \code{blackbt} object and prints a summary.
}

\usage{
  \method{summary}{blackbt}(object, \ldots)             
}

\arguments{
  \item{object}{ a \code{blackbt} output object. }
  \item{...}{ further arguments to \code{print}. }
}

\value{
    A summary of a \code{blackbt} object. For each dimension, reports all
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
'\link{blackbox_transpose}', '\link{LC1980}', '\link{plot.blackbt}', '\link{plotcdf.blackbt}', '\link{LC1980_bbt}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  LCdat <- LC1980[,-1]	# Dump the column of self-placements

  \donttest{ 
  LC1980_bbt <- blackbox_transpose(LCdat, missing=c(0,8,9), dims=3, 
    minscale=5, verbose=TRUE)
  }
  ### 'LC1980_bbt' can be retrieved quickly with: 
  data(LC1980_bbt)

  summary(LC1980_bbt)
  plot(LC1980_bbt)
}

\keyword{ multivariate }
