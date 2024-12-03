\name{plot.cdf}
\alias{plot.cdf}
\title{ Aldrich-McKelvey Coordinate Cumulative Distribution Plot }
\description{ 
  \code{plot.aldmck} reads an \code{aldmck} object and plots the cumulative distribution
  of the respondents and stimuli.
}

\usage{
\method{plot}{cdf}(x, align=NULL, xlim=c(-2,2), \ldots)
}

\arguments{
  \item{x}{ an \code{aldmck} output object. }
  \item{align}{ integer, the x-axis location that stimuli names should be aligned to
	If set to NULL, it will attempt to guess a location. }
  \item{xlim}{ vector of length 2, fed to the \code{plot} function as the \code{xlim}
	argument, which sets the minimum and maximum range of the x-axis.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A plot of the empirical cumulative distribution of the respondent ideal points,
    along with the locations of the stimuli. If no self-placements were specified
    during estimation, no graphical plots will appear.
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
'\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}', '\link{plot.aldmck}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), verbose=TRUE)
  
  summary(result)
  
  plot.cdf(result)
}

\keyword{ multivariate }
