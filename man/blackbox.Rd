\name{blackbox}
\alias{blackbox}
\title{ Blackbox Scaling }
\description{ 
  \code{blackbox} is a function that takes a matrix of survey data in which individuals
	place themselves on continuous scales across multiple issues, and locates those
	citizens in a spatial model of voting. Mathematically, this function generalizes
	the singular value of a matrix to cases in which there is missing data in the
	matrix. Scales generated using perceptual data (i.e. scales of legislator locations
	using liberal-conservative rankings by survey respondents) should instead use
	the \code{blackbox_transpose} function in this package instead.
}

\usage{
  blackbox(data,missing=NULL,verbose=FALSE,dims=1,minscale)
}

\arguments{
  \item{data}{ matrix of numeric values containing the issue scale data.  Respondents should
	be organized on rows, and stimuli on columns. It is helpful, though
	not necessary, to include row names and column names. }
  \item{missing}{ vector or matrix of numeric values, sets the missing values for the data.
	Observations with missing data are discarded before analysis.  If input is a
	vector, then the vector is assumed to contain the missing value codes for all
	the data.  If the input is a matrix, it must be of dimension p x q, where p
	is the maximum number of missing values and q is the number of columns in the
	data.  Each column of the inputted matrix then specifies the missing data values
	for the respective variables in data.  If null (default), no missing values are in
	the data.}
  \item{verbose}{ logical, indicates whether \code{aldmck} should print out detailed
	output when scaling the data.}
  \item{dims}{integer, specifies the number of dimensions to be estimated. }
  \item{minscale}{ integer, specifies the minimum number of responses a respondent needs
	needs to provide to be used in the scaling.} 
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


\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{rosentha@princeton.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}
}

\references{

Keith T. Poole (1998) ``Recovering a Basic Space From a Set of Issue Scales.'' American
Journal of Political Science. 42(3), 954-993.

}

\seealso{
  '\link{Issues1980}', '\link{summary.blackbox}'.
}

\examples{

### Loads issue scales from the 1980 NES.
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
result <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
summary(result)
}
\keyword{ multivariate }
