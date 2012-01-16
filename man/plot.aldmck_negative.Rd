\name{plot.aldmck_negative}
\alias{plot.aldmck_negative}
\title{ Aldrich-McKelvey Negative Coordinate Distribution Plot }
\description{ 
  \code{plot.aldmck_negative} reads an \code{aldmck} object and plots the probability distribution
  of the respondents and stimuli with negative weights.
}

\usage{
\method{plot}{aldmck_negative}(x, xlim=c(-2,2), \ldots)
}

\arguments{
  \item{x}{ an \code{aldmck} output object. }
  \item{xlim}{ vector of length 2, fed to the \code{plot} function as the \code{xlim}
	argument, which sets the minimum and maximum range of the x-axis.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A plot of the probability distribution of the respondent ideal points, along with
    the locations of the stimuli. If no negative weights exist, either because respondent
    self-placements are not specified, or because all weights are positive, a plot indicating
    this in text is given.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{rosentha@princeton.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}

}

\seealso{
  '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}', '\link{plot.cdf}', '\link{plot.aldmck}' 
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
summary(result)
plot.aldmck_negative(result,xlim=c(-1.5,1.5))
}
\keyword{ multivariate }
