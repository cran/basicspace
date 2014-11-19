\name{LC1980_bbt}
\alias{LC1980_bbt}

\title{
  Blackbox Transpose Estimate, 1980 Liberal-Conservative Scales.
}

\description{
  Blackbox-Transpose estimates from Liberal-Conservative 7-point scales from the 1980
  National Election Study. Estimates in 3 dimensions.
}

\usage{
data(LC1980_bbt)
}

\value{
    An object of class \code{blackbt}.

  \item{stimuli}{ vector of data frames of length dims. Each data frame presents results for
	estimates from that dimension (i.e. x$stimuli[[2]] presents results for dimension 2).
	Each row contains data on a separate stimulus, and each data frame includes the
	following variables:
    \itemize{
    \item{\code{N}}{Number of respondents who ranked this stimulus.}
    \item{\code{coord1D}}{Location of the stimulus in the first dimension. If viewing
	the results for a higher dimension, higher dimension results will appear as
	coord2D, coord3D, etc.}
    \item{\code{R2}}{The percent variance explained for the stimulus. This increases as
	more dimensions are estimated.}
    }                 
   }

  \item{individuals}{ vector of data frames of length dims. Each data frame presents results for
	estimates from that dimension (i.e. x$stimuli[[2]] presents results for dimension 2).
	Individuals that are discarded from analysis due to the minscale constraint are NA'd out.
	Each row contains data on a separate stimulus, and each data frame includes the
	following variables:
    \itemize{
    \item{\code{c}}{Estimate of the individual intercept.}
    \item{\code{w1}}{Estimate of the individual slope. If viewing the results for a higher
	dimension, higher dimension results will appear as w2, w3, etc.}
    \item{\code{R2}}{The percent variance explained for the respondent. This increases as
	more dimensions are estimated.}
    }                 
   }
  \item{fits}{ A data frame of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{SSE}}{Sum of squared errors.}
    \item{\code{SSE.explained}}{Explained sum of squared error.}
    \item{\code{percent}}{Percentage of total variance explained.}
    \item{\code{SE}}{Standard error of the estimate, with formula provided in the article cited below.}
    \item{\code{singular}}{Singluar value for the dimension.}
    }                 
  \item{Nrow}{ Number of rows/stimuli.}
  \item{Ncol}{ Number of columns used in estimation. This may differ from the data set due to
	columns discarded due to the minscale constraint.}
  \item{Ndata}{ Total number of data entries.}
  \item{Nmiss}{ Number of missing entries.}
  \item{SS_mean}{ Sum of squares grand mean.}
  \item{dims}{ Number of dimensions estimated.}
}

\source{

  American national Election Study. \url{http://www.electionstudies.org/}

}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{plotcdf.blackbt}', '\link{LC1980}', '\link{plot.blackbt}', '\link{summary.blackbt}', '\link{blackbox_transpose}'.
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

\keyword{datasets}
