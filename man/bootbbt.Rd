\name{bootbbt}
\alias{bootbbt}

\title{
  Blackbox Transpose Bootstrap of 1980 Liberal-Conservative Scales.
}

\description{
  Output from 10 bootstrap trials of LC1980 data. Included to allow the example to run
  sufficiently quickly to pass CRAN guidelines.
}

\usage{
  data(bootbbt)
}

\value{
  See '\link{boot_blackbt}'.
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
'\link{LC1980}', '\link{boot_blackbt}', '\link{plot.boot_blackbt}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  LCdat <- LC1980[,-1]	#Dump the column of self-placements

  \donttest{ 
    bootbbt <- boot_blackbt(LCdat, missing=c(0,8,9), dims=1, minscale=8, iter=10)
  }
  ### 'LC1980_bbt' can be retrieved quickly with: 
  data(bootbbt)

  plot.boot_blackbt(bootbbt)
}

\keyword{datasets}
