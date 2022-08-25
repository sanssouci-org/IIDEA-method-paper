# - - - - - - - - - - - - -
# packages and scripts
# - - - - - - - - - - - - -
library("sanssouci")
library("future.apply")
library("matrixStats")
library("cherry")
library("tibble")

source("scripts/utils/test_JER_control.R")
source("scripts/utils/add_signal.R")
source("scripts/utils/format_power.R")

# - - - - - - - - - - - - -
# data set
# - - - - - - - - - - - - -
technology <- "simulation"
n <- 60
m <- 10000
K <- 99 
rho <- c(0, 0.1, 0.4, 0.6, 0.8, 1)

ds_name <- "block-cov"

mu <- rep(0, m)
truthFp <- c(1, sort(sample(1:m, K, replace = FALSE)), m)
# truthS <- sort(sample(1:K, m, replace = T))
truthF <- rep(0, m)
for(i in 1:(length(truthFp)-1)){
  indF <- truthFp[i]:truthFp[i+1]
  truthF[indF] <- i
  # Sigma[indF, indF] <- Sigma[indF,indF] + rho
}

rowTestFUN <- sanssouci::rowPearsonCorrelationTests

# - - - - - - - - - - - - -
# parameters
# - - - - - - - - - - - - -
alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
B <- 1e3          # number of permutations for adaptive methods
nb_exp <- 1e3

snr <- 2
effect <- 1
nb_correlation = c(1, 2, 5, 10, 20, 50)

configs <- expand.grid(nb_correlation = nb_correlation)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_%s_correlation", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

for (r in rho){
  Sigma <- matrix(0, nrow = m, ncol = m)
  for (i in unique(truthF)){
    # indS <- which(truthS == i)
    indF <- which(truthF == i)
    Sigma[indF, indF] <- Sigma[indF, indF] + r
  }
  diag(Sigma) <- 1
  # image(Sigma)
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  X <- t(X)
  
  for (cc in seq_configs) {
    nb_cor <- configs[cc, "nb_correlation"]
    
    simname <- sprintf("%s_m=%s_nb_cor=%s_rho=%s_nb-exp=%s",
                       ds_name, m, nb_cor, r, nb_exp)
    print(simname)
    cat(cc, "/", nrow(configs), ":", simname, "\n")
    filename <- sprintf("%s.rds", simname)
    pathname <- file.path(path, filename)
    
    t0 <- Sys.time()
    # res <- lapply(1:nb_exp, FUN = function(i) {
    res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
      select <- sample(unique(truthF), nb_cor)
      non_zero = which(truthF %in% select)
      truth = rep(0, m)
      truth[non_zero] = effect
      eps = rnorm(n = n)
      prod_temp = truth %*% X
      noise_mag <- sqrt((prod_temp) %*% t(prod_temp)) / (snr * sqrt(eps%*% eps))
      Y = prod_temp + noise_mag[[1]] * eps
      
      tests <- rowTestFUN(X, Y)
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
        H = 1:m)

      ## check JER control and estimate power
      res_i <- test_JER_control_cont_cov(
        Y = X, groups = Y, truth = truth, 
        rowTestFUN = rowTestFUN, B = B, 
        alpha = alphas, selections = selections,
        verbose = TRUE)
      len <- sapply(selections, FUN = length) %>% 
        data.frame() %>% 
        tibble::rownames_to_column("selection") %>% 
        rename('length' = '.')
      list(level = tibble(exp = i, res_i$level),
           power = tibble(exp = i, res_i$power) %>% 
             left_join(., len, by = "selection"))
    })
    level <- Reduce(rbind, lapply(res, "[[", "level"))
    power <- Reduce(rbind, lapply(res, "[[", "power"))
    res <- list(level = level, power = power)
    saveRDS(res, file = pathname)
    print(simname)
  }
}
