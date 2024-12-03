\name{boot_blackbt}
\alias{boot_blackbt}
\title{ Bootstrap of Blackbox Transpose Scaling }
\description{ 
  \code{boot_blackbt} is a function automates the non-parametric bootstrapping of \code{blackbox_transpose}.
        The original function takes a matrix of perceptual data, such as liberal-conservative
        rankings of various stimuli, and recovers the true location of those stimuli in a spatial
        model. The bootstrap simply applies this estimator across multiple resampled data sets
        and stores the results of each iteration in a matrix.  These results can be used to estimate
        uncertainty for various parameters of interest, and can be plotted using the \code{plot.boot_blackbt}
        function.
}

\usage{
  boot_blackbt(data, missing=NULL, dims=1, dim.extract=dims, minscale, 
    iter=100, verbose=FALSE)
}

\arguments{
  \item{data}{ matrix of numeric values, containing the perceptual data.  Respondents should
	be organized on rows, and stimuli on columns. It is helpful, though
	not necessary, to include row names and column names. }
  \item{missing}{ vector or matrix of numeric values, sets the missing values for the data.
        NA values are always treated as missing regardless of what is set here.
	Observations with missing data are discarded before analysis.  If input is a
	vector, then the vector is assumed to contain the missing value codes for all
	the data.  If the input is a matrix, it must be of dimension p x q, where p
	is the maximum number of missing values and q is the number of columns in the
	data.  Each column of the inputted matrix then specifies the missing data values
	for the respective variables in data.  If null (default), no missing values are in
	the data other than the standard NA value.}
  \item{dims}{integer, specifies the number of dimensions to be estimated. }
  \item{dim.extract}{integer, specifies which dimension to extract results for the bootstrap from. }
  \item{minscale}{integer, specifies the minimum number of responses a respondent needs
	needs to provide to be used in the scaling.} 
  \item{iter}{integer, number of iterations the bootstrap should run for.}
  \item{verbose}{logical, indicates whether the progress of \code{blackbox_transpose} (at each 10th iteration) should be printed to the screen.}

}

\value{
    An object of class \code{boot_blackbt}. This is simply a matrix of dimensions iter x number of
    stimuli.  Each row stores the estimated stimuli locations for each iteration.

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
'\link{blackbox_transpose}', '\link{plot.boot_blackbt}'.
}


\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  LCdat <- LC1980[,-1]	#Dump the column of self-placements

  \donttest{ 
    bootbbt <- boot_blackbt(LCdat, missing=c(0,8,9), dims=1, 
      minscale=8, iter=10, verbose=FALSE)
  }
  ### 'LC1980_bbt' can be retrieved quickly with: 
  data(bootbbt)

  plot.boot_blackbt(bootbbt)
}

\keyword{ multivariate }
