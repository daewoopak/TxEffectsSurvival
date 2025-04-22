\name{wr.test}
\alias{wr.test}
\alias{wr.test.default}
\alias{print.wr.test}
\title{
The main function of the package provides various confidence intervals and testing procedures with event-specific win ratios
}
\description{
The function \code{wr.test} provides several confidence interval and testing procedures with the event-specific win ratios that are obtained on the terminal and non-terminal events. The following procedures are provided:
\itemize{
  \item Tests of the global null - testing the null hypothesis of no treatment effect on either the terminal event or the non-terminal event. A set of three tests are provided: the maximum test, the linear combination test, and the chi-squared test.
  \item Test of proportional hazards - testing the null hypothesis of the proportionality assumptions for the terminal event and the non-terminal event.
  \item Test of equal hazard ratios - testing the null hypothesis of equal hazard ratios for the terminal event and the non-terminal event when they both have proportional hazards.
  \item Confidence intervals of the non-terminal and terminal events respectively
  \item Confidence intervals of linear combinations of the non-terminal and terminal events, with either pre-determined or data-driven weights
}
The full details for these procedures are available in Yang et al. (2021).
}
\usage{
\method{wr.test}{default}(yh, hcen, yd, dcen, z, lin = c(0.5, 0.5), alpha = 0.05, repnum = 1E6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{for S4 method only.}
  \item{yh}{A numeric vector for time to the non-terminal event or censoring}
  \item{hcen}{Censoring indicator for the non-terminal event (event = 1, censored = 0)}
  \item{yd}{A numeric vector for time to the terminal event or censoring}
  \item{dcen}{Censoring indicator for the terminal event  (event = 1, censored = 0)}
  \item{z}{A numeric vector for the group indicator (treatment = 1, control = 0)}
  \item{lin}{A numeric vector of length 2 for the linear combination of the event-specific win ratios. The components must be non-negative values and added up to one. The first component is for the non-terminal event and the second is for the terminal event. The default is (0.5, 0.5).}
  \item{alpha}{The Significance level being used for confidence intervals. The default value is 0.05.}
  \item{repnum}{The number of replications for simulating bivariate normal distributions to obtain critical values corresponding to the alpha. The default value is 1E6.}
}


\value{
A S3 \code{wr.test} class object, which is a \code{list} with the following components:
\item{wr1}{The event specific win ratio for the non-terminal event}
\item{wr2}{The event specific win ratio for the terminal event}
\item{ci1t}{The confidence interval for the event specific win ratio for the non-terminal event}
\item{ci2t}{The confidence interval for the event specific win ratio for the terminal event}
\item{mxot}{The test statistic for the maximum test}
\item{pvalmxt}{The p-value for the maximum test}
\item{chi}{The test statistic for the chi-squre test}
\item{pvachi}{The p-value for the chi-squre test}
\item{lin}{The inputted vector for the linear combination of the event-specific win ratios}
\item{zvalin0}{The test statistic for the linear combination test}
\item{plin0}{The p-value for the linear combination test}
\item{wrlin0}{The weighted average win ratio with the inputted vector \code{lin}}
\item{cilin0}{The confidence interval for the weighted average win ratio with the inputted vector \code{lin}}
\item{lin_ar}{The data-driven linear combination}
\item{zvalint}{The test statistic for the data-driven Linear combination test}
\item{plintr}{The p-value for the data-driven Linear combination test}
\item{wrlinl}{The weighted average win ratio with the data-driven combination}
\item{cilint}{The confidence interval for the weighted average win ratio with the data-driven combination}
\item{mxph}{The test statistic for the test of proportional hazards}
\item{pvalph}{The p-value for the test of proportional hazards}
\item{zvaephl}{The test statistic for the test of equal hazard ratios}
\item{pvaephl}{The p-value for the test of equal hazard ratios}
}

\note{
Linear combination tests can be used to detect an overall effect, which is measured by using a weighted average of the win ratios of the terminal and non-terminal events, as considered in Yang and Troendle (2021). The weights can be either a data-driven weights or pre-determined weights. The pre-determined weights can be supplied with the \code{lin} argument.
}

\references{
Yang, S., Troendle, J., Pak, D., & Leifer, E. (2022). Event‐specific win ratios for inference with terminal and non‐terminal events. Statistics in medicine, 41(7), 1225-1241.

Yang, S., & Troendle, J. (2021). Event-specific win ratios and testing with terminal and non-terminal events. Clinical Trials, 18(2), 180-187.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
