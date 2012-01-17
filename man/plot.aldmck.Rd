\name{plot.aldmck}
\alias{plot.aldmck}
\title{ Aldrich-McKelvey Coordinate Distribution Plot }
\description{ 
  \code{plot.aldmck} reads an \code{aldmck} object and plots the probability distribution
  of the respondents and stimuli.
}

\usage{
\method{plot}{aldmck}(x, \ldots)
}

\arguments{
  \item{x}{ an \code{aldmck} output object. }
  \item{...}{ Other arguments to \code{plot}. }
}

\value{
    A plot of the probability distribution of the respondent ideal points, along with
    the locations of the stimuli. If no self-placements were specified during estimation,
    no graphical plots will appear.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{rosentha@princeton.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}

}

\seealso{
  '\link{aldmck}', '\link{LC1980}', '\link{summary.aldmck}', '\link{plot.AM}', '\link{plot.cdf}' '\link{plot.aldmck_negative}','\link{plot.aldmck_positive}'.
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
summary(result)
plot(result)
}
\keyword{ multivariate }
