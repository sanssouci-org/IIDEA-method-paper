# - - - - - - - - - - - - -
# packages and scripts
# - - - - - - - - - - - - -
library("sanssouci")
library("matrixStats")
library("future.apply")
library("tidyr")
library("dplyr")

source("scripts/utils/get_perm_lambda_FDP.R")
source("scripts/utils/format_power.R")

# - - - - - - - - - - - - -
# data set
# - - - - - - - - - - - - -
technology <- "RNAseq"
ds_name <- "BLCA"

data("RNAseq_blca", package = "sanssouci.data")
X <- RNAseq_blca
groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
rm(RNAseq_blca)

# filter out unexpressed genes
BLCA0 <- X / colSums(X) * 1e6
ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
if (length(ww) != 0) {
  X <- X[-ww, ]
}

# X <- head(X, 1234)

m <- nrow(X)
rowTestFUN <- sanssouci::rowWilcoxonTests

# - - - - - - - - - - - - -
# parameters
# - - - - - - - - - - - - -
alphas <- seq(from = 0, to = 1, by = 0.05) # target JER level
Bs <- c(100, 200, 500, 1000, 2000, 5000) # number of permutations
nb_exp <- 1000 # number of experiments

configs <- expand.grid(B = Bs)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_%s_permutation", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

# - - - - - - - - - - - - -
# experiments
# - - - - - - - - - - - - -
for (cc in seq_configs) {
  B <- configs[cc, "B"]

  simname <- sprintf(
    "%s_m=%s_B=%s_nb-exp=%s",
    ds_name, m,
    B, nb_exp
  )
  print(simname)
  cat(cc, "/", nrow(configs), ":", simname, "\n")
  filename <- sprintf("%s.rds", simname)
  pathname <- file.path(path, filename)

  t0 <- Sys.time()
  # res <- lapply(1:nb_exp, FUN = function(i) {
  res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
    print(i)
    tests <- rowTestFUN(X, groups)
    p_values <- tests$p.value
    m <- length(p_values)
    rk <- rank(p_values)
    selections <- list(
      # first_1 = which(rk %in% c(1:1)),
      first_10 = which(rk %in% c(1:10)),
      first_100 = which(rk %in% c(1:100)),
      first_1000 = which(rk %in% c(1:1000)),
      # first_5000 = which(rk %in% c(1:5000)),
      # first_10000 = which(rk %in% c(1:10000)),
      # first_15000 = which(rk %in% c(1:15000)),
      # BH_10 = which(p.adjust(p_values, method = "BH") < 0.10),
      BH_05 = which(p.adjust(p_values, method = "BH") < 0.05),
      p_05 = which(p_values < 0.05),
      # p_01 = which(p_values < 0.01),
      H = 1:m
    )
    ## check JER control and estimate power
    res_i <- get_perm_lambda_FDP(
      Y = X, groups = groups,
      rowTestFUN = rowTestFUN, B = B,
      alpha = alphas, selections = selections,
      verbose = TRUE
    )
    gc()
    list(
      lambda = tibble(exp = i, res_i$lambda),
      FDP = tibble(exp = i, res_i$FDP)
    )
  })
  lambda <- Reduce(rbind, lapply(res, "[[", "lambda"))
  FDP <- Reduce(rbind, lapply(res, "[[", "FDP"))
  res <- list(lambda = lambda, FDP = FDP)
  saveRDS(res, file = pathname)
}
