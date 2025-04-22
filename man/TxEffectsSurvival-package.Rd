\name{EventWinRatios-package}
\alias{EventWinRatios-package}
\docType{package}
\title{
  Event-Specific Win Ratios for Terminal and Non-Terminal Events
}
\description{
This package provides several confidence interval and testing procedures using event-specific win ratios for semi-competing risks data with non-terminal and terminal events, as developed in Yang et al. (2021). It allows users to utilize more data for improved inference procedures with semi-competing risks data. The event-specific win ratios were introduced in Yang and Troendle (2021), and the detailed definitions and asymptotic results for event-specific win ratios are available in Yang et al. (2021).

}
\details{
\tabular{ll}{
Package: \tab EventWinRatios\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2021-07-09\cr
License: \tab GPL (>= 3)\cr
}
}

\author{
Daewoo Pak and Song Yang

Maintainer: Daewoo Pak <heavyrain.pak@gmail.com>
}
\references{
Yang, S., Troendle, J., Pak, D., & Leifer, E. (2022). Event‐specific win ratios for inference with terminal and non‐terminal events. Statistics in medicine, 41(7), 1225-1241.

Yang, S., & Troendle, J. (2021). Event-specific win ratios and testing with terminal and non-terminal events. Clinical Trials, 18(2), 180-187.
}
\keyword{ package }
\seealso{
  \code{\link{wr.test}}
}
\examples{
library(EventWinRatios)
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
