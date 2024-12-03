\name{colombia}
\alias{colombia}

\title{
  2004 PELA Liberal-Conservative Scales.
}

\description{
  Liberal-Conservative 10-point scales from the University of Salamanca's Parliamentary Elites of Latin America (PELA) survey. Stored as a matrix of integers. The number 99 is a missing value. These data come from Sebastian Saiegh and are used in the paper and book cited below.
}

\usage{
data(colombia)
}

\value{
  The data is formatted as an integer matrix with the following elements.

  \item{colombia}{ matrix, containing reported placements of various stimuli on a
	10 point Liberal-Conservative scale:
    \describe{
    \item{\code{id}}{ Respondent ID.}
    \item{\code{party}}{ Respondent party.}
    \item{\code{departam}}{ Respondent district.}
    \item{\code{entrey}}{ Interviewer ID.}
    \item{\code{pl_uribista}}{ Placement of ``Partido Liberal Uribista'' on 10 point scale.}
    \item{\code{pl_oficial}}{ Placement of ``Partido Liberal Oficial'' on 10 point scale.}
    \item{\code{conservator}}{ Placement of ``Partido Conservador'' on 10 point scale.}
    \item{\code{polo}}{ Placement of ``Polo'' on 10 point scale.}
    \item{\code{union_cristiana}}{ Placement of ``Union Cristiana'' on 10 point scale.}
    \item{\code{salvation}}{ Placement of ``Salvacion'' on 10 point scale.}
    \item{\code{urine}}{ Placement of Mr. Uribe on 10 point scale.}
    \item{\code{antanas}}{ Placement of Mr. Antanas on 10 point scale.}
    \item{\code{gomez}}{ Placement of Mr. Gomez on 10 point scale.}
    \item{\code{garzon}}{ Placement of Garzon on 10 point scale.}
    \item{\code{holgin}}{ Placement of Holguin on 10 point scale.}
    \item{\code{rivera}}{ Placement of Rivera on 10 point scale.}
    \item{\code{self}}{ Respondent self placement on 10 point scale.}

    }
  }
}

\source{
Sebastian Saiegh. 2009. ``Recovering a Basic Space from Elite Surveys: Evidence from Latin America.'' \emph{Legislative Studies Quarterly} 34(1): 117-145. doi: 10.3162/036298009787500349

Sebastian Saiegh. 2011. \emph{Ruling By Statute: How Uncertainty and Vote-Buying Shape Lawmaking}. New York: Cambridge University Press. doi: 10.1017/CBO9780511842276
}

\references{
John H. Aldrich and Richard D. McKelvey. 1977. ``A Method of Scaling with Applications to the 1968 and 1972 Presidential Elections.'' \emph{American Political Science Review} 71(1): 111-130. doi: 10.2307/1956957

David A. Armstrong II, Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2021. \emph{Analyzing Spatial Models of Choice and Judgment}. 2nd ed. Statistics in the Social and Behavioral Sciences Series. Boca Raton, FL: Chapman & Hall/CRC. doi: 10.1201/9781315197609

Thomas R. Palfrey and Keith T. Poole. 1987. ``The Relationship between Information, Ideology, and Voting Behavior.'' \emph{American Journal of Political Science} 31(3): 511-530. doi: 10.2307/2111281

Keith T. Poole, Jeffrey B. Lewis, Howard Rosenthal, James Lo, and Royce Carroll. 2016. ``Recovering a Basic Space from Issue Scales in R.'' \emph{Journal of Statistical Software} 69(7): 1-21. doi:10.18637/jss.v069.i07

Keith T. Poole. 1998. ``Recovering a Basic Space From a Set of Issue Scales.'' \emph{American Journal of Political Science} 42(3): 954-993. doi: 10.2307/2991737
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}

Christopher Hare \email{cdhare@ucdavis.edu}
}

\seealso{
'\link{aldmck}', '\link{summary.aldmck}', '\link{plot.aldmck}', '\link{plot.cdf}'.
}

\examples{
  ### Loads the Liberal-Conservative scales from the 2004 PELA survey.
  data(colombia)
  tmp <- colombia[,c(5:8,12:17)]

  result <- aldmck(data=tmp, polarity=7, respondent=10, missing=c(99), verbose=TRUE)
  summary(result)
  plot.cdf(result)
}

\keyword{datasets}
