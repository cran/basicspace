\name{summary.blackbt}
\alias{summary.blackbt}
\title{ Blackbox-Transpose Summary }
\description{ 
  \code{summary.blackbt} reads an \code{blackbt} object and prints a summary.
}

\usage{
  summary.blackbt(object, ...)             
}

\arguments{
  \item{object}{ a \code{blackbt} output object. }
  \item{...}{ further arguments to \code{print}. }
}

\value{
    A summary of a \code{blackbt} object. For each dimension, reports all
    stimuli with coordinates, individuals used for scaling, and fit. Also
    summarizes number of rows, columns, total data entries, number of missing
    entries, percent missing data, and sum of squares.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{rosentha@princeton.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{blackbox_transpose}', '\link{LC1980}', '\link{plot.blackbt}', '\link{plotcdf.blackbt}'.
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
LCdat=LC1980[,-1]	#Dump the column of self-placements
result <- blackbox_transpose(LCdat,missing=c(0,8,9),dims=3,minscale=5,verbose=TRUE)
plot(result)
par(ask=TRUE)
plotcdf.blackbt(result)
summary(result)
}
\keyword{ multivariate }
