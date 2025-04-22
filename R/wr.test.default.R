wr.test.default <- function(yh, hcen, yd, dcen, z, lin = c(0.5, 0.5), alpha = 0.05, repnum = 1E6, tau_r = 0.9, ...) {

  trp <- max(c(yh, yd)) + 0.001
  trpn <- trp * tau_r

  match.call()

  yh <- as.numeric(yh)
  hcen <- as.numeric(hcen)
  yd <- as.numeric(yd)
  dcen <- as.numeric(dcen)
  z <- as.numeric(z)
  lin_raw <- lin
  z_alpha <- qnorm(1 - alpha/2)
  n <- length(z)

  set.seed(0)
  norr <- matrix(rnorm(repnum * 2), ncol = 2)
  ak <- norr[,1]

  ## data summary
  sum_table <- fn_ds(hcen, dcen, z)

  ## non-terminal
  gwrsim_h <- fn_wr(yh, hcen, z, ak, trp)
  wrall1 <- gwrsim_h$wrall
  zvaall1 <- gwrsim_h$zvaall
  stall1 <- gwrsim_h$stall
  zvaph1 <- gwrsim_h$zvaph
  zvata1 <- gwrsim_h$zvat4
  stnd1 <- gwrsim_h$stnd
  tqt1 <- gwrsim_h$tqt
  tdif1 <- gwrsim_h$tdif
  mhat01 <- gwrsim_h$mhat0
  mhatp101 <- gwrsim_h$mhatp10
  mhatp201 <- gwrsim_h$mhatp20
  mhatph01 <- gwrsim_h$mhatph

  ## terminal
  gwrsim_d <- fn_wr(yd, dcen, z, ak, trp)
  wrall2 <- gwrsim_d$wrall
  zvaall2 <- gwrsim_d$zvaall
  stall2 <- gwrsim_d$stall
  zvaph2 <- gwrsim_d$zvaph
  zvata2 <- gwrsim_d$zvat4
  stnd2 <- gwrsim_d$stnd
  tqt2 <- gwrsim_d$tqt
  tdif2 <- gwrsim_d$tdif
  mhat02 <- gwrsim_d$mhat0
  mhatp102 <- gwrsim_d$mhatp10
  mhatp202 <- gwrsim_d$mhatp20
  mhatph02 <- gwrsim_d$mhatph

  # CI for wr
  wr1 <- wrall1[1]
  ci1 <- c(wr1 - z_alpha * stnd1, wr1 + z_alpha * stnd1)

  lef1 <- tqt1 - z_alpha * stnd1 * tdif1;
  rit1 <- tqt1 + z_alpha * stnd1 * tdif1;
  ci1t <- c(log(2/(1-exp(-exp(lef1)))-1), log(2/(1-exp(-exp(rit1)))-1));

  wr2 <- wrall2[1]
  ci2 <- c(wr2 - z_alpha * stnd2, wr2 + z_alpha * stnd2)

  lef2 <- tqt2 - z_alpha * stnd2 * tdif2;
  rit2 <- tqt2 + z_alpha * stnd2 * tdif2;
  ci2t <- c(log(2/(1-exp(-exp(lef2)))-1), log(2/(1-exp(-exp(rit2)))-1));

  # CI for all
  mhatall <- cbind(mhat01, mhat02);
  mhatp1all <- cbind(mhatp101, mhatp102);
  mhatp2all <- cbind(mhatp201, mhatp202);
  wralla <- rbind(wrall1, wrall2);
  stalla <- rbind(stall1, stall2);

  ### p-values for wr
  rho <- sum(mhat01 * mhat02)/prod(sqrt(colSums(cbind(mhat01, mhat02)^2)));
  cmat <- matrix(c(1, rho, rho, 1), 2, 2)
  svd_cmat <- svd(cmat)
  rsig <- svd_cmat$u %*% diag(sqrt(svd_cmat$d)) %*% t(svd_cmat$v)

  sa <- norr %*% rsig
  # maximum test (21)
  mxo <- max(abs(c(zvaall1[1], zvaall2[1])))
  pvalmx <- mean(pmax(abs(sa[,1]), abs(sa[,2])) > mxo)

  # maximum test (21) - transformed one
  mxot <- max(abs(c(zvata1[1], zvata2[1])))
  pvalmxt <- mean(pmax(abs(sa[,1]), abs(sa[,2])) > mxot)

  mxpval <- c(pvalmx, pvalmxt)

  # Chisquare test (24) - transformed one
  crn <- 1 + 15 / max(c(n, 100))
  vecr <- c(zvata1, zvata2)/crn
  chi <- as.numeric(t(vecr) %*% solve(cmat) %*% vecr)
  pvachi <- pchisq(chi, df = 2, lower.tail = FALSE)
  # pvachi <- mean(rowSums(norr^2) > chi)

  # p-values for maxph (26)
  rhoph <- sum(mhatph01 * mhatph02)/prod(sqrt(colSums(cbind(mhatph01, mhatph02)^2)))
  zph12 <- c(zvaph1, zvaph2)
  cmath <- matrix(c(1, rhoph, rhoph, 1), 2, 2)

  svd_cmath <- svd(cmath)
  rsigh <- svd_cmath$u %*% diag(sqrt(svd_cmath$d)) %*% t(svd_cmath$v)
  sah <- norr %*% rsigh;
  mxph <- max(abs(c(zvaph1, zvaph2)));
  pvalph <- mean(pmax(abs(sah[,1]), abs(sah[,2])) > mxph)

  # linear
  ctv <- 1;
  # hrh <- 1; hr <- 1
  # wlm1 <- hrh
  # wlm2 <- hr
  # wmlm <- c(wlm1, wlm2)

  glinsim_l0 <- fn_glinsim(lin, wralla, mhatall, mhatp1all, mhatp2all, ak, ctv, n, z_alpha)
  wrlin0 <- glinsim_l0$wrlin[1]
  linci <- glinsim_l0$linci
  zvalin <- glinsim_l0$zvalin
  stlin <- glinsim_l0$stlin
  cilin0 <- linci[1,]
  plin0 <- mean(abs(ak) > abs(zvalin[1]))
  zas <- c(zvalin[1], stlin)

  ar <- stall1[1]^2 + stall2[1]^2 - 2*stall1[1]*stall2[1]*rho;
  ar <- (stall2[1]^2 - stall1[1]*stall2[1]*rho)/ar;

  if (ar < 0) {
    ar <- 0
  } else {
    if (ar > 1) {
      ar <- 1
    }
  }

  lin_ar <- c(ar, 1-ar);
  ctv <- 1;
  glinsim <- fn_glinsim(lin_ar, wralla, mhatall, mhatp1all, mhatp2all, ak, ctv, n, z_alpha)
  wrlinl <- glinsim$wrlin[1]
  lincil <- glinsim$lincil
  zvalinl <- glinsim$zvalinl
  plintr <- mean(abs(ak) > abs(zvalinl[1]))
  cilint <- lincil[1, ]

  # equal hazard test
  gephsim <- fn_gephsim(lin = c(-1, 1), wralla, mhatall, n, ak)
  zvaephl <- gephsim$zvaephl
  pvaephl <- gephsim$pvaephl

  #
  # covp <- as.numeric((linci[1,1] < lin %*% c(wlm1, wlm2)) * (linci[1,2] > lin %*% c(wlm1, wlm2)))
  # tlint0 <- (abs(zvalin[1]) > z_alpha)
  # cilin0 <- c(linci[1,1], linci[1,2])
  # len0 <- cilin0[2] - max(cilin0[1], 0)
  #
  #
  #
  #
  # covpt <- as.numeric((lincil[1,1] < lin_ar %*% c(wlm1, wlm2)) * (lincil[1,2] > lin_ar %*% c(wlm1, wlm2)))
  #
  # tlint <- (abs(zvalinl[1]) > z_alpha);
  # cilint <- c(lincil[1,1], lincil[1,2]);
  # lent <- cilint[2] - max(cilint[1], 0);
  #

  ### NEW approaches
  newfit <- fn_tntinf(yh, hcen, yd, dcen, zo, lin, norr, trp, trpn, ak)

  y_first <- pmin(yh, yd)
  d_first <- dcen + (1 - dcen) * hcen
  z_first <- zo
  logrank <- fn_logrank(y_first, d_first, z_first)
  wrtest0 <- fn_wr0(yh, hcen, yd, dcen, zo)

  results <- list()
  results$alpha <- alpha
  results$sum_table <- sum_table
  results$wr1 <- wr1 #Non-terminal event: win.ratio
  results$wr2 <- wr2
  results$ci1t <- ci1t #Non-terminal event: win.ratio confidence interval
  results$ci2t <- ci2t #terminal event: win.ratio confidence interval

  results$mxot <- mxot # test stat for Tests of the Global Null Hypothesis
  results$pvalmxt <- pvalmxt # p-value for Tests of the Global Null Hypothesis
  results$chi <- chi # chisquare test stat
  results$pvachi <- pvachi # chisquare test p-value

  results$lin <- lin_raw # input lin
  results$zvalin0 <- zvalin[1] #Linear combination test*: test stat
  results$plin0 <- plin0 #Linear combination test*: pvalue

  results$wrlin0 <- wrlin0 #Linear combination*: average win.ratio
  results$cilin0 <- cilin0 #Linear combination*: average win.ratio CI

  results$zvalint <- zvalinl[1] #Data-driven Linear combination test*: test.stat
  results$plintr <- plintr # Data-driven Linear combination test*: p-value

  results$wrlinl <- wrlinl # Data-driven Linear combination*: average win.ratio
  results$cilint <- cilint # Data-driven Linear combination*: average win.ratio CI
  results$lin_ar <- lin_ar # Data-driven Linear combination*: average win.ratio linear

  results$mxph <- mxph # Test of Proportional Hazards test stat
  results$pvalph <- pvalph #Test of Proportional Hazards pvalue
  results$zvaephl <- zvaephl #Test of Equal Hazard Ratios test stat
  results$pvaephl <- pvaephl # Test of Equal Hazard Ratios pvalue

  results$newfit <- newfit # new methods
  results$logrank <- logrank #logrank
  results$wrtest0 <- wrtest0

  class(results) <- "wr.test"

  results

}

