\name{LC1980}
\alias{LC1980}

\title{
  1980 Liberal-Conservative Scales.
}

\description{
  Liberal-Conservative 7-point scales from the 1980 National Election Study. Includes
  (in order) self-placement, and rankings of Carter, Reagan, Kennedy, Anderson,
  Republican party, Democratic Party. Stored as a matrix of integers. The numbers 0,
  8, and 9 are considered to be missing values.
}

\usage{
data(LC1980)
}

\value{

  The data is formatted as an integer matrix with the following elements.

  \item{LC1980}{ matrix, containing reported placements of various stimuli on a
	7 point Liberal-Conservative scale:
    \itemize{
    \item{\code{Self}}{ Self-placement on 7 point scale.}
    \item{\code{Carter}}{ Placement of Carter on 7 point scale.}
    \item{\code{Reagan}}{ Placement of Reagan on 7 point scale.}
    \item{\code{Kennedy}}{ Placement of Kennedy on 7 point scale.}
    \item{\code{Anderson}}{ Placement of Anderson on 7 point scale.}
    \item{\code{Republicans}}{ Placement of Republicans on 7 point scale.}
    \item{\code{Democrats}}{ Placement of Democrats on 7 point scale.}
    }
  }
}

\source{

  American national Election Study. \url{http://www.electionstudies.org/}

  Also availble from Keith Poole. \url{ http://voteview.edu/aldmck.htm }
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{lo@uni-mannheim.de}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@ucla.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{aldmck}', '\link{summary.aldmck}', '\link{plot.aldmck}', '\link{plot.cdf}'.
}

\examples{
### Loads and scales the Liberal-Conservative scales from the 1980 NES.
data(LC1980)
result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9),verbose=TRUE)
summary(result)
par(ask=TRUE)
plot.AM(result,xlim=c(-1.5,1.5))
plot.cdf(result)
}

\keyword{datasets}
