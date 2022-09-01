test_JER_control <- function(Y, groups, truth, rowTestFUN = sanssouci::rowWelchTests,
                             B = 100, alpha = 0.05,
                             selections = NULL, verbose = FALSE) {
  obj_i <- SansSouci(Y = Y, groups = groups)

  # Oracle predictions: depend on experiment (via gene order) but not on alpha!
  obj_oracle <- obj_i
  obj_oracle$input$truth <- truth
  obj_oracle <- fit(obj_oracle, rowTestFUN = rowTestFUN, family = "Oracle", alpha = NA)
  FP_oracle <- predict(obj_oracle, what = "FP", all = TRUE)$bound

  TP_oracle <- sapply(selections, FUN = function(sel) {
    predict(obj_oracle, S = sel, what = "TP", all = FALSE)
  })

  level0 <- power0 <- NULL ## Simes (no calibration)
  level <- power <- NULL ## Simes + single-step calibration
  level_sd <- power_sd <- NULL ## Simes + step-down calibration
  level_cherry <- power_cherry <- NULL ## cherry/ARI

  cal0 <- fit(obj_i,
    B = 0, rowTestFUN = rowTestFUN,
    alpha = alpha[1],
    family = "Simes"
  )

  cal <- fit(obj_i,
    B = B, rowTestFUN = rowTestFUN,
    alpha = alpha[1],
    family = "Simes", max_steps_down = 0
  )
  piv_stat <- cal$output$piv_stat # does not depend on alpha[1] because single-step
  m <- nHyp(obj_i)

  p_values <- pValues(cal0)
  hom <- hommelFast(p_values)

  for (aa in seq(along = alpha)) {
    ## Simes + calibration (single-step)
    lambda <- stats::quantile(piv_stat, alpha[aa], type = 1)
    thr <- t_linear(lambda, 1:m, nHyp(cal)) ## hack to avoid recalc. pivotal stat
    cal$output$thr <- thr

    FP <- predict(cal, what = "FP", all = TRUE)$bound
    valid_bound <- all(FP >= FP_oracle)
    level[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)

    TP <- sapply(selections, FUN = function(sel) {
      predict(cal, S = sel, what = "TP", all = FALSE)
    })
    pow <- format_power(TP / TP_oracle)
    power[[aa]] <- tibble(alpha = alpha[aa], pow)

    ## Simes + calibration (step-down)
    cal_sd <- fit(cal,
      B = B, rowTestFUN = rowTestFUN,
      alpha = alpha[aa], family = "Simes"
    )

    FP <- predict(cal_sd, what = "FP", all = TRUE)$bound
    valid_bound <- all(FP >= FP_oracle)
    level_sd[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)

    TP_sd <- sapply(selections, FUN = function(sel) {
      predict(cal, S = sel, what = "TP", all = FALSE)
    })
    pow_sd <- format_power(TP_sd / TP_oracle)
    power_sd[[aa]] <- tibble(alpha = alpha[aa], pow_sd)

    ## Simes + no calibration
    thr <- t_linear(alpha[aa], 1:m, nHyp(cal)) ## hack to avoid recalc. pivotal stat
    cal0$output$thr <- thr

    FP0 <- predict(cal0, what = "FP", all = TRUE)$bound
    valid_bound <- all(FP0 >= FP_oracle)
    level0[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)

    TP0 <- sapply(selections, FUN = function(sel) {
      predict(cal0, S = sel, what = "TP", all = FALSE)
    })
    pow0 <- format_power(TP0 / TP_oracle)
    power0[[aa]] <- tibble(alpha = alpha[aa], pow0)

    ## cherry/ARI
    TP_cherry <- curveSimes(hommel = hom, alpha = alpha[aa], plot = FALSE)
    FP_cherry <- 1:m - TP_cherry
    valid_bound <- all(FP_cherry >= FP_oracle)
    level_cherry[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)

    TP_cherry <- sapply(selections, FUN = function(sel) {
      pickSimes(hommel = hom, select = sel, alpha = alpha[aa], silent = TRUE)
    })
    pow_cherry <- format_power(TP_cherry / TP_oracle)
    power_cherry[[aa]] <- tibble(alpha = alpha[aa], pow_cherry)
  }

  level <- tibble(method = "Simes + single-step calibration", Reduce(rbind, level))
  level0 <- tibble(method = "Simes (parametric)", Reduce(rbind, level0))
  level_sd <- tibble(method = "Simes + step-down calibration", Reduce(rbind, level_sd))
  level_cherry <- tibble(method = "cherry/ARI", Reduce(rbind, level_cherry))

  power <- tibble(method = "Simes + single-step calibration", Reduce(rbind, power))
  power0 <- tibble(method = "Simes (parametric)", Reduce(rbind, power0))
  power_sd <- tibble(method = "Simes + step-down calibration", Reduce(rbind, power_sd))
  power_cherry <- tibble(method = "cherry/ARI", Reduce(rbind, power_cherry))

  list(
    level = rbind(level, level0, level_sd, level_cherry),
    power = rbind(power, power0, power_sd, power_cherry)
  )
}
