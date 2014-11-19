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

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}

}

\seealso{
  '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}', '\link{plot.aldmck}'.
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)

summary(result)
par(ask=TRUE)
plot.AM(result,xlim=c(-1.5,1.5))
plot.cdf(result)
}
\keyword{ multivariate }
