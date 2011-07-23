\name{plot.aldmck_positive}
\alias{plot.aldmck_positive}
\title{ Aldrich-McKelvey Positive Coordinate Distribution Plot }
\description{ 
  \code{plot.aldmck_positive} reads an \code{aldmck} object and plots the probability distribution
  of the respondents and stimuli with positive weights.
}

\usage{
plot.aldmck_positive(x, xlim=c(-2,2), ...)
}

\arguments{
  \item{x}{ an \code{aldmck} output object. }
  \item{xlim}{ vector of length 2, fed to the \code{plot} function as the \code{xlim}
	argument, which sets the minimum and maximum range of the x-axis.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A plot of the probability distribution of the respondent ideal points, along with
    the locations of the stimuli. If no weights exist because respondent self-placements
    are not specified, a plot indicating this in text is given.
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
plot.aldmck_positive(result,xlim=c(-1.5,1.5))
}
\keyword{ multivariate }
