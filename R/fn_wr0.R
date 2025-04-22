fn_wr0 <- function(yh, hcen, yd, dcen, zo) {

  n <- length(zo)

  w1 <- l1 <- w2 <- l2 <- w1o <- l1o <- sig2o <- sig1o <- sig1 <- c()
  for (i in 1:n) {
    w2ij <- dcen * (yd[i] >= yd)
    l2ij <- dcen[i] * (yd[i] <= yd)

    w1ij <- hcen * (yh[i] >= yh)
    l1ij <- hcen[i] * (yh[i] <= yh)

    w1ijo <- (1 - w2ij) * (1 - l2ij) * w1ij
    l1ijo <- (1 - w2ij) * (1 - l2ij) * l1ij

    w1[i] <- sum(zo[i] * (1 - zo) * w1ij)
    l1[i] <- sum(zo[i] * (1 - zo) * l1ij)

    w2[i] <- sum(zo[i] * (1 - zo) * w2ij)
    l2[i] <- sum(zo[i] * (1 - zo) * l2ij)

    w1o[i] <- sum(zo[i] * (1 - zo) * w1ijo)
    l1o[i] <- sum(zo[i] * (1 - zo) * l1ijo)

    sig2o[i] <- mean((zo[i] - zo) * (w2ij - l2ij))
    sig1o[i] <- mean((zo[i] - zo) * (1 - w2ij) * (1 - l2ij) * (w1ij - l1ij))
    sig1[i] <- mean((zo[i] - zo) * (w1ij - l1ij))
  }

  sigd <- sqrt(mean((sig1o + sig2o)^2) / n)
  wd <- (-sum(l2 + l1o) + sum(w2 + w1o)) / n^2
  zwr <- wd / sigd

  # pvawr <- mean(abs(ak) > abs(zwr))
  pvawr <- 2 * pnorm(abs(zwr), lower.tail = FALSE)

  results <- list()
  results$stat_wr <- zwr
  results$pval_wr <- pvawr

  results


}

