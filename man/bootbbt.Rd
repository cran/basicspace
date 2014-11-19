\name{bootbbt}
\alias{bootbbt}

\title{
  Blackbox Transpose Bootstrap of 1980 Liberal-Conservative Scales.
}

\description{
  Output from 10 bootstrap trials of LC1980 data. Included to allow the example to run
  sufficiently quickly to pass CRAN guidelines.
}

\usage{
data(bootbbt)
}

\value{

  See '\link{boot_blackbt}'.

}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{LC1980}', '\link{boot_blackbt}', '\link{plot.boot_blackbt}'.
}

\examples{


data(LC1980)

data=LC1980[,-1]


# Not run to save time, but loaded object is the output
 # bootbbt <- boot_blackbt(data, missing=c(0,8,9), dims=1, minscale=8, iter=10)
 data("bootbbt")

 plot.boot_blackbt(bootbbt)


}

\keyword{datasets}
