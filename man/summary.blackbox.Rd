\name{summary.blackbox}
\alias{summary.blackbox}
\title{ Blackbox Summary }
\description{ 
  \code{summary.blackbox} reads an \code{blackbox} object and prints a summary.
}

\usage{
  \method{summary}{blackbox}(object, \ldots)             
}

\arguments{
  \item{object}{ a \code{blackbox} output object. }
  \item{...}{ further arguments to \code{print}. }
}

\value{
    A summary of a \code{blackbox} object. For each dimension, reports all
    stimuli with coordinates, individuals used for scaling, and fit. Also
    summarizes number of rows, columns, total data entries, number of missing
    entries, percent missing data, and sum of squares.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{blackbox}', '\link{Issues1980}'
}

\examples{
### Loads issue scales from the 1980 NES.
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

### This command conducts estimates, which we instead load using data()
# Issues1980_bb <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
data(Issues1980_bb)

summary(Issues1980_bb)
}
\keyword{ multivariate }
