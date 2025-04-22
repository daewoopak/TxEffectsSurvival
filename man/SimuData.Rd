\name{SimuData}
\alias{SimuData}
\docType{data}
\title{A simulated semi-competing risks data set with non-terminal events and terminal events}
\description{
\code{SimuData} is a simulated data set that is used to show the utility of this package. Under the semi-competing risks data setup, each subject has a non-terminal and a terminal event that may be right censored. Furthermore, the terminal event may censor the non-terminal event.
}
\usage{data("SimuData")}
\format{
  A data frame with 3500 observations on the following 5 variables.
  \describe{
    \item{\code{yh}}{Time to the non-terminal event of interest or censoring (by an independent censoring variable or a terminal event)}
    \item{\code{hcen}}{Censoring indicator for the non-teminal event (event = 1, censored = 0)}
    \item{\code{yd}}{Time to the terminal event of interest}
    \item{\code{dcen}}{Censoring indicator for the teminal event (event = 1, censored = 0)}
    \item{\code{z}}{Group indicator (treatment = 1, control = 0)}
  }
}

\examples{
data(SimuData)

# non-terminal events
yh <- SimuData$yh
hcen <- SimuData$hcen

# terminal events
yd <- SimuData$yd
dcen <- SimuData$dcen

# group indicator
z <- SimuData$z

# Win Ratio tests
result <- wr.test(yh, hcen, yd, dcen, z)
print(result)

}
\seealso{\code{\link{wr.test}}}
\keyword{datasets}
