\name{Issues1980_bb}
\alias{Issues1980_bb}

\title{
  Blackbox Estimate, 1980 NES Issue Scales.
}

\description{
  Blackbox estimates from issues scales from the 1980 National Election Study.
}

\usage{
data(Issues1980_bb)
}

\value{
    An object of class \code{blackbox}.

  \item{stimuli}{ vector of data frames of length dims. Each data frame presents results for
	estimates from that dimension (i.e. x$stimuli[[2]] presents results for dimension 2).
	Each row contains data on a separate stimulus, and each data frame includes the
	following variables:
    \itemize{
    \item{\code{N}}{Number of respondents who provided a response to this stimulus.}
    \item{\code{c}}{Stimulus intercept.}
    \item{\code{w1}}{Estimate of the stimulus weight on the first dimension. If viewing the
	results for a higher dimension, higher dimension results will appear as w2, w3, etc.}
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
    \item{\code{c1}}{Estimate of the individual intercept on the first dimension. If viewing the
	results for a higher dimension, higher dimension results will appear as c2, c3, etc.}
    }                 
   }
  \item{fits}{ A data frame of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{SSE}}{Sum of squared errors.}
    \item{\code{SSE.explained}}{Explained sum of squared error.}
    \item{\code{percent}}{Percentage of total variance explained.}
    \item{\code{SE}}{Standard error of the estimate, with formula provided on pg. 973 of the article cited below.}
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

  American National Election Study. \url{http://www.electionstudies.org/}

}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{Issues1980}', '\link{summary.blackbox}', '\link{plot.blackbox}'.
}

\examples{

### Loads issue scales from the 1980 NES.
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

### This command conducts estimates, which we instead load using data()
#Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
data(Issues1980_bb)

summary(Issues1980_bb)

}

\keyword{datasets}
