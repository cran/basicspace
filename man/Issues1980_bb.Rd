\name{Issues1980_bb}
\alias{Issues1980_bb}

\title{
  Blackbox Estimate, 1980 ANES Issue Scales.
}

\description{
  Blackbox estimates from issues scales from the 1980 American National Election Study.
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
    \describe{
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
    \describe{
    \item{\code{c1}}{Estimate of the individual intercept on the first dimension. If viewing the
	results for a higher dimension, higher dimension results will appear as c2, c3, etc.}
    }                 
   }
  \item{fits}{ A data frame of fit results, with elements listed as follows:}
    \describe{
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
  American National Election Studies (\url{https://electionstudies.org/})
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
'\link{Issues1980}', '\link{summary.blackbox}', '\link{plot.blackbox}'.
}

\examples{
  ### Loads issue scales from the 1980 ANES.
  data(Issues1980)
  Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
  Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

  \donttest{ 
  Issues1980_bb <- blackbox(Issues1980, missing=c(0,8,9), verbose=FALSE, 
    dims=3, minscale=8)
  }
  ### 'Issues1980_bb' can be retrieved quickly with: 
  data(Issues1980_bb)

  summary(Issues1980_bb)
}

\keyword{datasets}
