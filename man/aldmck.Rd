\name{aldmck}
\alias{aldmck}
\title{Aldrich-McKelvey Scaling}
\description{ 
  \code{aldmck} is a function that takes a matrix of perceptual data, such as
	liberal-conservative rankings of various stimuli, and recovers the true
	location of those stimuli in a spatial model. It differs from procedures
	such as \code{wnominate}, which instead use preference data to estimate
	candidate and citizen positions. The procedure here, developed by John
	Aldrich and Richard McKelvey in 1977, is restricted to estimating data
	with no missing values and only in one dimension.  Please refer to the
	\code{blackbox} and \code{blackbox_transpose} functions in this package for
	procedures that accomodate missing data and multidimensionality estimates.
}

\usage{
  aldmck(data, respondent=0, missing=NULL, polarity, verbose=FALSE)
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
  \item{verbose}{ logical, indicates whether \code{aldmck} should print out detailed
	output when scaling the data.}
}

\value{
    An object of class \code{aldmck}.

  \item{legislators}{ vector, containing the recovered locations of the stimuli. The names of
	the stimuli are attached if provided as column names in the argument \code{data},
	otherwise they are generated sequentiall as 'stim1', 'stim2', etc. }

  \item{respondents}{ matrix, containing the information estimated for each respondent.
	Observations which were discarded in the estimation for missing data purposes
	have been NA'd out:
    \describe{
    \item{\code{intercept}}{Intercept of perceptual distortion for respondent.}
    \item{\code{weight}}{Weight of perceptual distortion for respondent.}
    \item{\code{idealpt}}{Estimated location of the respondent. Note that these positions
	are still calculated for individuals with negative weights, so these may need
	to be discarded. Note that this will not be calculated if self-placements are
	not provided in the data.}
    \item{\code{selfplace}}{The self-reported location of the individual, copied from the
	\code{data} argument if \code{respondent} is not set to 0. }
    \item{\code{polinfo}}{Estimated political information of respondent, calculated as
	the correlation between the true and reported stimulus locations.  The validation
	of this measure is provided in the article by Palfrey and Poole in the references.
	Note that this measure is included even for respondents that were not used in the
	estimation. Individuals with negative weights have also been assigned a political
	information score of 0, rather than negative scores.}
    }                 
  }
  \item{eigenvalues}{ A vector of the eigenvalues from the estimation.}
  \item{AMfit}{ Ratio of overall variance to perceptions in scaled data divided by average
	variance. This measure of fit, described by Aldrich and McKelvey, measures the
	amount of reduction of the variance of the scaled over unscaled data.}
  \item{N}{ Number of respondents used in the estimation (i.e. had no missing data)}
  \item{N.neg}{ Number of cases with negative weights. Only calculated if respondent
	self-placements are specified, will equal 0 if not.}
  \item{N.pos}{ Number of cases with positive weights. Only calculated if respondent
	self-placements are specified, will equal 0 if not.}

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
	
	result <- aldmck(data=LC1980, polarity=2, respondent=1, 
		missing=c(0,8,9), verbose=FALSE)
	
	summary(result)
	plot.aldmck(result)
}

\keyword{ multivariate }
