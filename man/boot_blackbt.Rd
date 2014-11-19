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
  boot_blackbt(data, missing=NULL, dims=1, dim.extract=dims, minscale, iter=100)
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
  \item{minscale}{ integer, specifies the minimum number of responses a respondent needs
	needs to provide to be used in the scaling.} 
  \item{iter}{integer, number of iterations the bootstrap should run for.}
}

\value{
    An object of class \code{boot_blackbt}. This is simply a matrix of dimensions iter x number of
    stimuli.  Each row stores the estimated stimuli locations for each iteration.

}


\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{blackbox_transpose}', '\link{plot.boot_blackbt}'.
}

\examples{

data(LC1980)

data=LC1980[,-1]

 # Not run to save time, but loaded object is the output

# bootbbt <- boot_blackbt(data, missing=c(0,8,9), dims=1, minscale=8, iter=10)

data("bootbbt")


plot.boot_blackbt(bootbbt)


}
\keyword{ multivariate }
