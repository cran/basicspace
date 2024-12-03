\name{plot.boot_blackbt}
\alias{plot.boot_blackbt}
\title{ Bootstrapped Blackbox Transpose Stimulus Plots }
\description{ 
  \code{plot.boot_blackbt} reads an \code{boot_blackbt} object and plots a dotchart of the stimuli with estimated confidence intervals.
}



\usage{

\method{plot}{boot_blackbt}(x, \ldots)

}



\arguments{

  \item{x}{ an \code{boot_blackbt} output object. }

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
'\link{blackbox_transpose}', '\link{boot_blackbt}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  data <- LC1980[,-1]	#Dump the column of self-placements

  \donttest{ 
  bootbbt <- boot_blackbt(data, missing=c(0,8,9), dims=1, 
    minscale=8, iter=10)
  }
  ### 'bootbbt' can be retrieved quickly with: 
  data(bootbbt)

  plot.boot_blackbt(bootbbt)
}

\keyword{ multivariate }

