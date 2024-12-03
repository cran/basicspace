\name{Issues1980}
\alias{Issues1980}

\title{
  1980 Issues Scakes
}

\description{
  Issue scales from the 1980 National Election Study. The numbers 0, 8, and 9 are considered
  to be missing values, except for the two abortion scales, where '7' is also a missing value.
  Hence, it must be recoded as in the example shown below before scaling. The data is used
  as an example for blackbox().
}

\usage{
data(LC1980)
}

\value{

  The data is formatted as an numeric matrix with the following elements.

  \item{Issues}{ matrix, containing reported self-placements along various stimuli on a
	7 point Liberal-Conservative scales (with the exception of abortion scales, which
	are 4 point):
    \describe{
    \item{\code{libcon1}}{ Liberal-conservative self-placement on 7 point scale.}
    \item{\code{defense}}{ Defense spending self-placement on 7 point scale.}
    \item{\code{govserv}}{ Government service on 7 point scale.}
    \item{\code{inflation}}{ Importance of inflation self-placement on 7 point scale.}
    \item{\code{abortion1}}{ Attitude on abortion 4 point scale.}
    \item{\code{taxcut}}{ Support for tax cut on 7 point scale.}
    \item{\code{libcon2}}{ Liberal-conservative self-placement on 7 point scale.}
    \item{\code{govhelpmin}}{ Government aid on 7 point scale.}
    \item{\code{russia}}{ Attitude towards Russia on 7 point scale.}
    \item{\code{womenrole}}{ Role of women on 7 point scale.}
    \item{\code{govjobs}}{ Placement of Democrats on 7 point scale.}
    \item{\code{equalrights}}{ Support for equal rights on 7 point scale.}
    \item{\code{busing}}{ Opinion on busing on 7 point scale.}
    \item{\code{abortion2}}{ Another attitude on abortion on 4 point scale.}

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
David A. Armstrong II, Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2021. \emph{Analyzing Spatial Models of Choice and Judgment}. 2nd ed. Statistics in the Social and Behavioral Sciences Series. Boca Raton, FL: Chapman & Hall/CRC. doi: 10.1201/9781315197609

Keith T. Poole, Jeffrey B. Lewis, Howard Rosenthal, James Lo, and Royce Carroll. 2016. ``Recovering a Basic Space from Issue Scales in R.'' \emph{Journal of Statistical Software} 69(7): 1-21. doi:10.18637/jss.v069.i07

Keith T. Poole. 1998. ``Recovering a Basic Space From a Set of Issue Scales.'' \emph{American Journal of Political Science} 42(3): 954-993. doi: 10.2307/2991737
}

\seealso{
'\link{blackbox}', '\link{summary.blackbox}'.
}

\examples{
  ### Loads issue scales from the 1980 ANES.
  data(Issues1980)
  Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
  Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode

  \donttest{ 
  Issues1980_bb <- blackbox(Issues1980, missing=c(0,8,9), verbose=FALSE, 
    dims=3, minscale=8)
  }
  ### 'Issues1980_bb' can be retrieved quickly with: 
  data(Issues1980_bb)

  summary(Issues1980_bb)
}

\keyword{datasets}
