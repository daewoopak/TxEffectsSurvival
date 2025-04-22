fn_gephsim <- function(lin, wralla, mhatall, n, ak) {
  
  wrnon <- wralla[1,1]
  wrter <- wralla[2,1]
  
  qtn <- 1 - exp(-wrnon^2)
  qtt <- 1 - exp(-wrter^2)
  dwrl <- log(-log(qtn)) - log(-log(qtt))
  dtn <- 2 * wrnon * exp(-wrnon^2)/qtn/log(qtn)
  dtt <- 2 * wrter * exp(-wrter^2)/qtt/log(qtt)
  dtm <- diag(c(dtn, dtt))

  stl <- sqrt(mean((mhatall %*% dtm %*% lin)^2))/sqrt(n)

  zvaephl <- dwrl/stl
  pvaephl <- mean(abs(ak) > abs(zvaephl))

  res <- list()
  res$zvaephl <- zvaephl
  res$pvaephl <- pvaephl
  
  res
  
}