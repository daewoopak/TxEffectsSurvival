logit_ci <- function(tq, st, tdif) {
  lef <- tq - 1.96 * st * tdif
  rit <- tq + 1.96 * st * tdif
  log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)
}
