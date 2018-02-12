\name{plot.boot_aldmck}
\alias{plot.boot_aldmck}
\title{ Bootstrapped Aldrich-McKelvey Stimulus Plots }
\description{ 
  \code{plot.boot_aldmck} reads an \code{boot_aldmck} object and plots a dotchart of the stimuli with estimated confidence intervals.
}



\usage{

\method{plot}{boot_aldmck}(x, \ldots)

}



\arguments{

  \item{x}{ an \code{boot_aldmck} output object. }

  \item{...}{ other arguments to \code{plot}.}

}



\value{

    A dotchart of estimated stimulus locations, with 95 percent confidence intervals. Point estimates are estimates from the original data set.
}



\author{ 

Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}

}



\seealso{

   '\link{aldmck}', '\link{boot_aldmck}'.
}



\examples{


data(LC1980)
result <- boot_aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), iter=30)
plot(result)

}

\keyword{ multivariate }

