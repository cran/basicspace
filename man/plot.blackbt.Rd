\name{plot.blackbt}
\alias{plot.blackbt}
\title{ Blackbox Transpose Coordinate Distribution Plot }
\description{ 
  \code{plot.blackbt} reads an \code{blackbt} object and plots the probability distribution
  of the respondents and stimuli.
}

\usage{
plot.blackbt(x, xlim=c(-1,1), ...)
}

\arguments{
  \item{x}{ an \code{blackbt} output object. }
  \item{xlim}{ vector of length 2, fed to the \code{plot} function as the \code{xlim}
	argument, which sets the minimum and maximum range of the x-axis.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A plot of the probability distribution of the respondent ideal points, along with
    the locations of the stimuli.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{rosentha@princeton.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}

}

\seealso{
  '\link{blackbox_transpose}', '\link{LC1980}', '\link{plotcdf.blackbt}', '\link{summary.blackbt}'.
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
LCdat=LC1980[,-1]	#Dump the column of self-placements
result <- blackbox_transpose(LCdat,missing=c(0,8,9),dims=3,minscale=5,verbose=TRUE)
plot(result)
par(ask=TRUE)
plotcdf.blackbt(result)
summary(result)
}
\keyword{ multivariate }
