\name{plot.boot_blackbt}
\alias{plot.boot_blackbt}
\title{ Bootstrapped Blackbox Transpose Stimulus Plots }
\description{ 
  \code{plot.boot_blackbt} reads an \code{boot_blackbt} object and plots a dotchart of the stimuli with estimated confidence intervals.
}



\usage{

\method{plot}{boot_blackbt}(x, \ldots)

}



\arguments{

  \item{x}{ an \code{boot_blackbt} output object. }

  \item{...}{ other arguments to \code{plot}.}

}



\value{

    A dotchart of estimated stimulus locations, with 95 percent confidence intervals. Point estimates are estimates from the original data set.
}



\author{ 

Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}

}



\seealso{

   '\link{blackbox_transpose}', '\link{boot_blackbt}'.
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

