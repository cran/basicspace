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
    \describe{
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
    \describe{
    \item{\code{c}}{Estimate of the individual intercept.}
    \item{\code{w1}}{Estimate of the individual slope. If viewing the results for a higher
	dimension, higher dimension results will appear as w2, w3, etc.}
    \item{\code{R2}}{The percent variance explained for the respondent. This increases as
	more dimensions are estimated.}
    }                 
   }
  \item{fits}{ A data frame of fit results, with elements listed as follows:}
    \describe{
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
'\link{plotcdf.blackbt}', '\link{LC1980}', '\link{plot.blackbt}', '\link{summary.blackbt}', '\link{blackbox_transpose}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  LCdat <- LC1980[,-1]	#Dump the column of self-placements

  \donttest{ 
  LC1980_bbt <- blackbox_transpose(LCdat, missing=c(0,8,9), dims=3, 
    minscale=5, verbose=TRUE)
  }
  ### 'LC1980_bbt' can be retrieved quickly with: 
  data(LC1980_bbt)

  summary(LC1980_bbt)
  plot(LC1980_bbt)
}

\keyword{datasets}
