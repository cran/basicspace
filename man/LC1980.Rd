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
    \describe{
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
  American National Election Studies (\url{https://electionstudies.org/})
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
'\link{aldmck}', '\link{summary.aldmck}', '\link{plot.aldmck}', '\link{plot.cdf}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 1980 ANES.
  data(LC1980)
  
  result <- aldmck(data=LC1980, polarity=2, respondent=1, missing=c(0,8,9), verbose=TRUE)

  summary(result)
  plot(result)
}

\keyword{datasets}
