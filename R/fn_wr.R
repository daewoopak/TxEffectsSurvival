fn_wr <- function(yor, dor, zor, ak, trp) {

  n <- length(yor)
  o <- order(yor)
  oy <- yor[o]
  od <- dor[o]
  oz <- zor[o]

  zd1 <- which(oz == 0)
  zd2 <- which(oz == 1)
  n1 <- length(zd1)
  n2 <- length(zd2)
  z1 <- rep(1, n1)
  z2 <- rep(0, n2)
  y1 <- oy[zd1]
  y2 <- oy[zd2]
  d1 <- od[zd1]
  d2 <- od[zd2]
  yk1 <- n1:1
  yk2 <- n2:1
  ky <- n:1

  x <- exp(-cumsum(od/ky))
  xat1 <- x[zd1]
  xat2 <- x[zd2]

  yk2at1 <- rep(1, n1)
  for (ik in 1:n1) {
    yk2at1[ik] <- sum(y2 >= y1[ik])
  }

  yk1at2 <- rep(1, n2)
  for (ik in 1:n2) {
    yk1at2[ik] <- sum(y1 >= y2[ik])
  }

  trunk2 <- as.numeric(y2 <= trp)
  trunk1 <- as.numeric(y1 <= trp)
  # trunk2n <- as.numeric(y2 <= trpn)
  # trunk1n <- as.numeric(y1 <= trpn)

  f1 <- yk2at1 * trunk1 / n
  f2 <- yk1at2 * trunk2 / n

  los <- as.numeric(t(d2)%*%f2/n)
  win <- as.numeric(t(d1)%*%f1/n)
  wr <- los/win #

  ph11 <- n/(yk2at1 + yk1)
  ph12 <- n/(yk2 + yk1at2)
  ph21 <- xat1 * ph11
  ph22 <- xat2 * ph12

  losp1 <- as.numeric(t(d2) %*% (ph12 * f2)/n)
  winp1 <- as.numeric(t(d1) %*% (ph11 * f1)/n)
  wrp1 <- losp1/winp1

  losp2 <- as.numeric(t(d2) %*% (ph22 * f2)/n)
  winp2 <- as.numeric(t(d1) %*% (ph21 * f1)/n)
  wrp2 <- losp2/winp2

  wrall <- c(wr, wrp1, wrp2)

  ff1 <- -wr * f1
  ff2 <- f2

  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)
  stnd <- sqrt(mean(mhat^2))/win/sqrt(n)
  mhat0 <- mhat/win

  zva <- (wr - 1)/stnd

  qt <- (exp(wr)-1)/(exp(wr)+1)
  qt0 <- (exp(1)-1)/(exp(1)+1)
  tqt <- log(-log(qt))
  tdif <- 2*exp(wr)/log(qt)/qt/(1+exp(wr))^2
  zvat4 <- (tqt-log(-log(qt0)))/tdif/stnd
  pvalt <- mean(abs(ak) > abs(zvat4))
  zvapvalt <- c(zvat4, pvalt)

  ff1 <- -wrp1 * f1 * ph11
  ff2 <- f2 * ph12
  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)
  stndp1 <- sqrt(mean(mhat^2))/winp1/sqrt(n)
  mhatp10 <- mhat/winp1
  zvap1 <- (wrp1-1)/stndp1

  ff1 <- -wrp2 * f1 * ph21
  ff2 <- f2* ph22
  mhat <- fn_gmhat(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n)
  stndp2 <- sqrt(mean(mhat^2))/winp2/sqrt(n)
  mhatp20 <- mhat/winp2
  zvap2 <- (wrp2-1)/stndp2

  stall <- c(stnd, stndp1, stndp2)
  zvaall <- c(zva, zvap1, zvap2)

  mhatph <- mhatp20 - mhatp10
  stndph <- sqrt(mean(mhatph^2))/sqrt(n)
  zvaph <- (wrp2-wrp1)/stndph

  rk <- rep(0, n);
  for (i in 1:n) {
    rk[i] <- sum(yor[i] >= yor)
  }

  pval <- mean(abs(ak) > abs(zva))
  pvalp1 <- mean(abs(ak) > abs(zvap1))
  pvalp2 <- mean(abs(ak) > abs(zvap2))
  pvalph <- mean(abs(ak) > abs(zvaph))
  pvall <- c(pval, pvalp1, pvalp2)


  res <- list()
  res$wrall <- wrall
  res$zvaall <- zvaall
  res$stall <- stall
  res$zvaph <- zvaph
  res$zvat4 <- zvat4
  res$stnd <- stnd
  res$tqt <- tqt
  res$tdif <- tdif
  res$mhat0 <- mhat0[rk]
  res$mhatp10 <- mhatp10[rk]
  res$mhatp20 <- mhatp20[rk]
  res$mhatph <- mhatph[rk]
  res$pval <- pval
  res$pvalp1 <- pvalp1
  res$pvalp2 <- pvalp2
  res$pvalph <- pvalph
  res$pvall <- pvall

  res


}

