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
    \itemize{
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

  American national Election Study. \url{http://www.electionstudies.org/}

  Also availble from Keith Poole. \url{ http://voteview.edu/}
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@princeton.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '  '\link{blackbox}', '\link{summary.blackbox}'.
}

\examples{
### Loads issue scales from the 1980 NES.
data(Issues1980)
Issues1980[Issues1980[,"abortion1"]==7,"abortion1"] <- 8	#missing recode
Issues1980[Issues1980[,"abortion2"]==7,"abortion2"] <- 8	#missing recode
result <- blackbox(Issues1980,missing=c(0,8,9),verbose=FALSE,dims=3,minscale=8)
summary(result)
}

\keyword{datasets}
