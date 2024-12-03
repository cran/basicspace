\name{boot_aldmck}
\alias{boot_aldmck}
\title{ Bootstrap of Aldrich-McKelvey Scaling }
\description{ 
  \code{boot_aldmck} is a function automates the non-parametric bootstrapping of \code{aldmck}.
        The original function takes a matrix of perceptual data, such as liberal-conservative
        rankings of various stimuli, and recovers the true location of those stimuli in a spatial
        model. The bootstrap simply applies this estimator across multiple resampled data sets
        and stores the results of each iteration in a matrix.  These results can be used to estimate
        uncertainty for various parameters of interest, and can be plotted using the \code{plot.boot_aldmck}
        function.
}

\usage{
  boot_aldmck(data, respondent = 0, missing=NULL, polarity, iter=100)
}

\arguments{
  \item{data}{ matrix of numeric values, containing the perceptual data.  Respondents should
	be organized on rows, and stimuli on columns. It is helpful, though
	not necessary, to include row names and column names. }
  \item{respondent}{ integer, specifies the column in the data matrix of the stimuli that contains
	the respondent's self-placement on the scale. Setting respondent = 0 specifies
	that the self-placement data is not available.  Self-placement data is not
	required to estimate the locations of the stimuli, but is required if recovery
	of the respondent ideal points, or distortion parameters is desired. Note that
	no distortion parameters are estimated in AM without self-placements because they
	are not needed, see equation (24) in Aldrich and McKelvey (1977) for proof.}
  \item{missing}{ vector or matrix of numeric values, sets the missing values for the data.
        NA values are always treated as missing regardless of what is set here.
	Observations with missing data are discarded before analysis.  If input is a
	vector, then the vector is assumed to contain the missing value codes for all
	the data.  If the input is a matrix, it must be of dimension p x q, where p
	is the maximum number of missing values and q is the number of columns in the
	data.  Each column of the inputted matrix then specifies the missing data values
	for the respective variables in data. If null (default), no missing values are in
	the data other than the standard NA value.}
  \item{polarity}{integer, specifies the column in the data matrix of the stimuli that is
	to be set on the left side (generally this means a liberal) }
  \item{iter}{ integer, is the number of iterations the bootstrap should run for.}
}

\value{
    An object of class \code{boot_aldmck}. This is simply a matrix of dimensions iter x number of
    stimuli.  Each row stores the estimated stimuli locations for each iteration.

}


\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{rosentha@princeton.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}

Christopher Hare \email{cdhare@ucdavis.edu}
}

\references{
John H. Aldrich and Richard D. McKelvey. 1977. ``A Method of Scaling with Applications to the 1968 and 1972 Presidential Elections.'' \emph{American Political Science Review} 71(1): 111-130. doi: 10.2307/1956957

David A. Armstrong II, Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2021. \emph{Analyzing Spatial Models of Choice and Judgment}. 2nd ed. Statistics in the Social and Behavioral Sciences Series. Boca Raton, FL: Chapman & Hall/CRC. doi: 10.1201/9781315197609

Thomas R. Palfrey and Keith T. Poole. 1987. ``The Relationship between Information, Ideology, and Voting Behavior.'' \emph{American Journal of Political Science} 31(3): 511-530. doi: 10.2307/2111281

Keith T. Poole, Jeffrey B. Lewis, Howard Rosenthal, James Lo, and Royce Carroll. 2016. ``Recovering a Basic Space from Issue Scales in R.'' \emph{Journal of Statistical Software} 69(7): 1-21. doi:10.18637/jss.v069.i07

Keith T. Poole. 1998. ``Recovering a Basic Space From a Set of Issue Scales.'' \emph{American Journal of Political Science} 42(3): 954-993. doi: 10.2307/2991737
}

\seealso{
'\link{LC1980}', '\link{summary.aldmck}',  '\link{plot.aldmck}', '\link{plot.cdf}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  
  result <- boot_aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), iter=30)

  plot(result)
}

\keyword{ multivariate }
