get_perm_lambda_FDP <- function(Y, groups, rowTestFUN = sanssouci::rowWelchTests,
                                B = 100, alpha = 0.05,
                                selections = NULL, verbose = FALSE) {
  obj_i <- SansSouci(Y = Y, groups = groups)


  lambda <- FDP <- NULL
  tic <- Sys.time()
  ## Simes + single step calibration
  cal <- fit(obj_i,
    B = B, rowTestFUN = rowTestFUN,
    alpha = alpha[1],
    family = "Simes", max_steps_down = 0
  )
  toc <- Sys.time()
  m <- nHyp(obj_i)

  for (aa in seq(along = alpha)) {
    ## Simes + calibration (step-down)
    cal_sd <- fit(cal,
      B = B, rowTestFUN = rowTestFUN,
      alpha = alpha[aa], family = "Simes"
    )

    lambda_a <- cal_sd[["output"]][["lambda"]]
    lambda[[aa]] <- tibble(alpha = alpha[aa], lambda = lambda_a, time = toc - tic)

    FDPa <- sapply(selections, FUN = function(sel) {
      predict(cal_sd, S = sel, what = c("FDP"), all = FALSE)
    })
    FDPa <- format_power(FDPa, values_to = "FDP")
    FDP[[aa]] <- tibble(alpha = alpha[aa], FDPa, time = toc - tic)
  }
  lambda <- tibble(method = "Simes + step-down calibration", Reduce(rbind, lambda))

  FDP <- tibble(method = "Simes + step-down calibration", Reduce(rbind, FDP))

  list(
    lambda = rbind(lambda),
    FDP = rbind(FDP)
  )
}
