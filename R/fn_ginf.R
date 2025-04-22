# yor <- yh; dor <- hcen
# trp <- 0.9 * max(c(yh, yd))
fn_ginf <- function(yor, dor, zo, trp, trpn, ak) {

  I <- order(yor)
  oy <- yor[I]
  od <- dor[I]
  oz <- zo[I]

  n <- length(yor)
  ky <- n:1

  zd1 <- which(oz == 0)
  zd2 <- which(oz == 1)
  n1 <- length(zd1)
  n2 <- length(zd2)
  n <- n1 + n2

  z1 <- rep(0, n1)
  z2 <- rep(1, n2)

  y1 <- oy[zd1]
  y2 <- oy[zd2]
  d1 <- od[zd1]
  d2 <- od[zd2]

  yk1 <- n1:1
  yk2 <- n2:1

  la1 <- cumsum(d1 / yk1)
  la2 <- cumsum(d2 / yk2)

  yk2at1 <- sapply(y1, function(y) sum(y2 >= y))
  yk1at2 <- sapply(y2, function(y) sum(y1 >= y))

  trunk2 <- as.numeric(y2 <= trp)
  trunk1 <- as.numeric(y1 <= trp)
  trunk2n <- as.numeric(y2 <= trpn)
  trunk1n <- as.numeric(y1 <= trpn)

  f1 <- yk2at1 * trunk1 / n
  f2 <- yk1at2 * trunk2 / n

  los <- sum(d2 * f2) / n
  win <- sum(d1 * f1) / n
  wr <- los / win

  ph11 <- n / (yk2at1 + yk1)
  ph12 <- n / (yk2 + yk1at2)
  losp1 <- sum(d2 * (ph12 * f2)) / n
  winp1 <- sum(d1 * (ph11 * f1)) / n
  wrp1 <- losp1 / winp1

  ntau2 <- sum(y2 < trpn)
  dytau2 <- c(y2[2:ntau2], trpn) - y2[1:ntau2]
  losl <- sum(la2[1:ntau2] * dytau2)

  ntau1 <- sum(y1 < trpn)
  dytau1 <- c(y1[2:ntau1], trpn) - y1[1:ntau1]
  winl <- sum(la1[1:ntau1] * dytau1)
  wrl <- losl / winl

  itchc <- sum((1 - exp(-la1[1:ntau1])) * dytau1)
  itcht <- sum((1 - exp(-la2[1:ntau2])) * dytau2)
  nwr <- itcht / itchc

  wrall <- c(wr, wrp1, wrl, nwr)

  # z-value 계산
  ff1 <- -wr * f1
  ff2 <- f2

  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)

  stnd <- sqrt(mean(mhat^2)) / win / sqrt(n)
  mhat0 <- mhat / win
  zva <- (wr - 1) / stnd

  qt <- (exp(wr) - 1) / (exp(wr) + 1)
  qt0 <- (exp(1) - 1) / (exp(1) + 1)
  tqt <- log(-log(qt))
  tdif <- 2 * exp(wr) / log(qt) / qt / (1 + exp(wr))^2
  zvat <- (tqt - log(-log(qt0))) / tdif / stnd

  # wrp1
  ff1 <- -wrp1 * f1 * ph11
  ff2 <- f2 * ph12
  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)
  stndp1 <- sqrt(mean(mhat^2)) / winp1 / sqrt(n)
  mhatp10 <- mhat / winp1
  zvap1 <- (wrp1 - 1) / stndp1

  qtp1 <- (exp(wrp1) - 1) / (exp(wrp1) + 1)
  tqtp1 <- log(-log(qtp1))
  tdifp1 <- 2 * exp(wrp1) / log(qtp1) / qtp1 / (1 + exp(wrp1))^2
  zvap1t <- (tqtp1 - log(-log(qt0))) / tdifp1 / stndp1

  # wrl
  g1 <- n * (trpn - y1) / yk1 * trunk1n
  g2 <- n * (trpn - y2) / yk2 * trunk2n
  ff1 <- -wrl * g1
  ff2 <- g2
  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)
  stndl <- sqrt(mean(mhat^2)) / winl / sqrt(n)
  mhatg0 <- mhat / winl
  zval <- (wrl - 1) / stndl

  qtl <- (exp(wrl) - 1) / (exp(wrl) + 1)
  tqtl <- log(-log(qtl))
  tdifl <- 2 * exp(wrl) / log(qtl) / qtl / (1 + exp(wrl))^2
  zvalt <- (tqtl - log(-log(qt0))) / tdifl / stndl

  # nwr
  la1l <- c(0, la1[1:(n1-1)])
  phi1 <- (y1 - c(0, y1[1:(n1-1)])) * exp(-la1l)
  phi1 <- cumsum(phi1)
  phi1tau <- phi1[ntau1] + exp(-la1[ntau1]) * (trpn - y1[ntau1])
  tch1 <- n * (phi1tau - phi1) / yk1 * trunk1n

  la2l <- c(0, la2[1:(n2-1)])
  phi2 <- (y2 - c(0, y2[1:(n2-1)])) * exp(-la2l)
  phi2 <- cumsum(phi2)
  phi2tau <- phi2[ntau2] + exp(-la2[ntau2]) * (trpn - y2[ntau2])
  tch2 <- n * (phi2tau - phi2) / yk2 * trunk2n

  ff1 <- -nwr * tch1
  ff2 <- tch2
  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)
  stndrm <- sqrt(mean(mhat^2)) / itchc / sqrt(n)
  mhatrm0 <- mhat / itchc
  zvarm <- (nwr - 1) / stndrm

  qtn <- (exp(nwr) - 1) / (exp(nwr) + 1)
  tqtn <- log(-log(qtn))
  tdifn <- 2 * exp(nwr) / log(qtn) / qtn / (1 + exp(nwr))^2
  zvarmt <- (tqtn - log(-log(qt0))) / tdifn / stndrm

  stall <- c(stnd, stndp1, stndl, stndrm)
  zvaall <- c(zva, zvap1, zval, zvarm)

  # p-value (normal approximation)
  pval <- mean(abs(ak) > abs(zva))
  pvalp1 <- mean(abs(ak) > abs(zvap1))
  pvall0 <- mean(abs(ak) > abs(zval))
  pvalrm <- mean(abs(ak) > abs(zvarm))
  pvall <- c(pval, pvalp1, pvall0, pvalrm)

  results <- list()
  results$wrall <- wrall
  results$stall <- stall
  results$zva <- zva
  results$zvat <- zvat
  results$zvap1 <- zvap1
  results$zval <- zval
  results$zvarm <- zvarm
  results$zvat <- zvat
  results$zvap1t <- zvap1t
  results$zvalt <- zvalt
  results$zvarmt <- zvarmt
  results$wrl <- wrl
  results$stndl <- stndl
  results$nwr <- nwr
  results$stndrm <- stndrm
  results$tqt <- tqt
  results$stnd <- stnd
  results$tdif <- tdif
  results$tqtp1 <- tqtp1
  results$stndp1 <- stndp1
  results$tdifp1 <- tdifp1
  results$tqtl <- tqtl
  results$tdifl <- tdifl
  results$tqtn <- tqtn
  results$tdifn <- tdifn
  results$mhat0 <- mhat0
  results$mhatp10 <- mhatp10
  results$mhatg0 <- mhatg0
  results$mhatrm0 <- mhatrm0

  results

}


