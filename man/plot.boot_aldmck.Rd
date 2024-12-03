\name{plot.boot_aldmck}
\alias{plot.boot_aldmck}
\title{ Bootstrapped Aldrich-McKelvey Stimulus Plots }
\description{ 
  \code{plot.boot_aldmck} reads an \code{boot_aldmck} object and plots a dotchart of the stimuli with estimated confidence intervals.
}



\usage{

\method{plot}{boot_aldmck}(x, \ldots)

}



\arguments{

  \item{x}{ an \code{boot_aldmck} output object. }

  \item{...}{ other arguments to \code{plot}.}

}



\value{

    A dotchart of estimated stimulus locations, with 95 percent confidence intervals. Point estimates are estimates from the original data set.
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
'\link{aldmck}', '\link{boot_aldmck}'.
}



\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  result <- boot_aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), iter=30)
  
  plot(result)
}

\keyword{ multivariate }

