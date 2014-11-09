\name{plotcdf.blackbt}
\alias{plotcdf.blackbt}
\title{ Blackbox Transpose Coordinate Cumulative Distribution Plot }
\description{ 
  \code{plotcdf.blackbt} reads an \code{blackbt} object and plots the cumulative distribution
  of the respondents and stimuli.
}

\usage{
  plotcdf.blackbt(x, align=NULL, xlim=c(-1.2,1), \ldots)
}

\arguments{
  \item{x}{ an \code{blackbt} output object. }
  \item{align}{ integer, the x-axis location that stimuli names should be aligned to
	If set to NULL, it will attempt to guess a location. }
  \item{xlim}{ vector of length 2, fed to the \code{plot} function as the \code{xlim}
	argument, which sets the minimum and maximum range of the x-axis.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A plot of the empirical cumulative distribution of the respondent ideal points,
    along with the locations of the stimuli.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}

}

\seealso{
  '\link{blackbox_transpose}', '\link{LC1980}', '\link{plot.blackbt}', '\link{summary.blackbt}', '\link{LC1980_bbt}'.
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
LCdat=LC1980[,-1]	#Dump the column of self-placements

### This command conducts estimates, which we instead load using data()
#LC1980_bbt <- blackbox_transpose(LCdat,missing=c(0,8,9),dims=3,minscale=5,verbose=TRUE)
data(LC1980_bbt)

plot(LC1980_bbt)
par(ask=TRUE)
plotcdf.blackbt(LC1980_bbt)
summary(LC1980_bbt)
}
\keyword{ multivariate }
