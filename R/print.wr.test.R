print.wr.test <- function(x, ...) {

  ci_level <- round((1 - x$alpha) * 100, 2)
  digit <- paste("%.", max(3, getOption("digits") - 4), "f", sep = "")
  digitp <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")

  display_p <- function(sv) {
    ifelse(as.numeric(sv) < 0.0001, "< 0.0001", paste("= ", sv, sep = ""))
  }

  wr1 <- sprintf(digit, x$wr1)
  wr2 <- sprintf(digit, x$wr2)
  # ci1 <- sprintf(digit, x$ci1)
  # ci2 <- sprintf(digit, x$ci2)
  ci1t <- sprintf(digit, x$ci1t)
  ci2t <- sprintf(digit, x$ci2t)

  mxo <- sprintf(digit, x$mxo)
  pvalmx <- display_p(sprintf(digitp, x$pvalmx))
  mxot <- sprintf(digit, x$mxot)
  pvalmxt <- display_p(sprintf(digitp, x$pvalmxt))

  chi <- sprintf(digit, x$chi)
  pvachi <- display_p(sprintf(digitp, x$pvachi))

  lin <- round(x$lin, 2)
  wrlin0 <- sprintf(digit, x$wrlin0)
  zvalin0 <- sprintf(digit, x$zvalin0)
  cilin0 <- sprintf(digit, x$cilin0)
  plin0 <- display_p(sprintf(digitp, x$plin0))

  lin_ar <- round(x$lin_ar, 2)
  wrlinl <- sprintf(digit, x$wrlinl)
  zvalint <- sprintf(digit, x$zvalint)
  cilint <- sprintf(digit, x$cilint)
  plintr <- display_p(sprintf(digitp, x$plintr))

  mxph <- sprintf(digit, x$mxph)
  pvalph <- display_p(sprintf(digitp, x$pvalph))
  zvaephl <- sprintf(digit, x$zvaephl)
  pvaephl <- display_p(sprintf(digitp, x$pvaephl))

  ## newfit
  new_wralla <- x$newfit$wralla
  new_ci1 <- x$newfit$ci1
  new_ci2 <- x$newfit$ci2
  new_cil <- x$newfit$cil
  new_mx <- x$newfit$mx
  new_pvalmx <- x$newfit$pvalmx
  new_zvalin <- x$newfit$zvalin
  new_pvalinall <- x$newfit$pvalinall
  new_wrlin <- x$newfit$wrlin
  new_stlinall <- x$newfit$stlinall

  new_mxf <- sprintf(digit, new_mx)
  new_pvalmxf <- display_p(sprintf(digitp, new_pvalmx))
  new_wrlinf <- sprintf(digit, new_wrlin)
  new_pvalinallf <- display_p(sprintf(digitp, new_pvalinall))

  new_wrlin_ar <- x$newfit$wrlin_ar
  new_pwrlin_ar <- x$newfit$peswrlog
  new_wrlin_ar_f <- sprintf(digit, new_wrlin_ar)
  new_pwrlin_ar_f <- display_p(sprintf(digitp, new_pwrlin_ar))

  new_stat_wr0 <- x$wrtest0$stat_wr
  new_pval_wr0 <- x$wrtest0$pval_wr
  new_stat_wr0_f <- sprintf(digit, new_stat_wr0)
  new_pval_wr0_f <- display_p(sprintf(digitp, new_pval_wr0))


  logrank_stat <- x$logrank$stat
  logrank_pval <- x$logrank$pval

  logrank_statf <- sprintf(digit, logrank_stat^2)
  logrank_pvalf <- display_p(sprintf(digitp, logrank_pval))

  sum_table <- x$sum_table
  colnames(sum_table) <- c("Treatment", "Control", "All")
  rownames(sum_table) <- c("   H only", "   H & D", "   D only")

  len <- 81
  cat(strrep("=", len), "\n")
  cat("I. -------------- < Summary of Observed Events > -------------------------------\n")
  print(sum_table)
  cat(strrep("-", len), "\n\n")

  cat("II. ------------- < Tests of the Global Null Hypothesis > -----------------------\n")

  cat(sprintf("\n>> Linear combination test [l = (%.1f, %.1f)]\n", lin[1], lin[2]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ESWR", new_wrlinf[1], new_pvalinallf[1]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "LRGRE", new_wrlinf[2], new_pvalinallf[2]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ICH", new_wrlinf[3], new_pvalinallf[3]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ITCH", new_wrlinf[4], new_pvalinallf[4]))

  cat(sprintf("\n>> Linear combination test with data-driven weights* [l = (%.2f, %.2f)]\n", lin_ar[1], lin_ar[2]))
  cat("*(log transformation used)\n")
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ESWR", new_wrlinf[1], new_pvalinallf[1]))

  cat("\n>> Maximum test (log-log transformation used)\n")
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ESWR", new_mxf[1], new_pvalmxf[1]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "LRGRE", new_mxf[2], new_pvalmxf[2]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ICH", new_mxf[3], new_pvalmxf[3]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ITCH", new_mxf[4], new_pvalmxf[4]))

  cat(paste("\n>> Chi-squred test: statistic = ", chi, ", p-value ", pvachi, sep = ""), "\n")
  cat(paste(">> Win ratio test: Chisq.statistic = ", new_stat_wr0_f, ", p-value ", new_pval_wr0_f, sep = ""), "\n")
  cat(paste(">> Log-rank test using time to first event: statistic = ", logrank_statf, ", p-value ", logrank_pvalf, sep = ""), "\n")
  cat(strrep("-", len), "\n\n")

  cat(paste("III. ----------- < Confidence Intervals (", ci_level, "% CI*) > -----------------------------", sep = ""), "\n")
  cat("*(loglog transformation used)\n\n")
  cat(">> Non-terminal event\n")
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ESWR", new_wralla[1, 1], new_ci1[3, 1], new_ci1[3, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "LRGRE", new_wralla[1, 2], new_ci1[4, 1], new_ci1[4, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ICH", new_wralla[1, 3], new_ci1[5, 1], new_ci1[5, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ITCH", new_wralla[1, 4], new_ci1[6, 1], new_ci1[6, 2]))

  cat("\n>> Terminal event\n")
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ESWR", new_wralla[2, 1], new_ci2[3, 1], new_ci2[3, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "LRGRE", new_wralla[2, 2], new_ci2[4, 1], new_ci2[4, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ICH", new_wralla[2, 3], new_ci2[5, 1], new_ci2[5, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ITCH", new_wralla[2, 4], new_ci2[6, 1], new_ci2[6, 2]))

  cat(sprintf("\n>> Linear combinations [l = (%.1f, %.1f)]\n", lin[1], lin[2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ESWR", new_wrlin[1], new_cil[1, 1], new_cil[1, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "LRGRE", new_wrlin[2], new_cil[2, 1], new_cil[2, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ICH", new_wrlin[3], new_cil[3, 1], new_cil[3, 2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ITCH", new_wrlin[4], new_cil[4, 1], new_cil[4, 2]))

  cat(sprintf("\n>> Data-driven combinations [l = (%.2f, %.2f)]\n", lin_ar[1], lin_ar[2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ESWR",
              as.numeric(wrlinl), as.numeric(cilint[1]), as.numeric(cilint[2])))
  cat(strrep("-", len), "\n\n")

   cat("IV. ------------ < Test of Proportional Hazards > -------------------------------\n")
  cat(paste(">> test.stat = ", mxph, ", p-value ", pvalph, sep = ""), "\n")
  cat(strrep("-", len), "\n\n")

  cat("V. ------------- < Test of Equal Hazard Ratios > --------------------------------\n")
  cat(paste(">> test.stat = ", zvaephl, ", p-value ", pvaephl, sep = ""), "\n")
  cat(strrep("-", len), "\n")
  cat(strrep("=", len))
}
