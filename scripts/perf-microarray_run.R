# - - - - - - - - - - - - -
# packages and scripts
# - - - - - - - - - - - - -
library("sanssouci")
library("future.apply")
library("matrixStats")
library("cherry")
library("tibble")
library("GSEABenchmarkeR")

source("scripts/utils/test_JER_control.R")
source("scripts/utils/add_signal.R")
source("scripts/utils/format_power.R")

# - - - - - - - - - - - - -
# data set
# - - - - - - - - - - - - -
technology <- "microarray"
rowTestFUN <- sanssouci::rowWelchTests
ds_name <- "GSE19188" # should be available from geo2kegg (via GSEABenchmarkeR)

source("scripts/utils/load_microarray_data.R") # loads 'X0' and 'groups'
m <- nrow(X0)
data_set <- sprintf("%s_m=%s", ds_name, m)
print(data_set)

path <- sprintf("results/diff-expr_%s_perf", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

# - - - - - - - - - - - - -
# parameters
# - - - - - - - - - - - - -
alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
B <- 1000          # number of permutations for adaptive methods
nb_exp <- 1000     # number of experiments

pi0s <- c(0.5, 0.8, 0.9, 1)
SNR <- c(0, 1, 2, 5)
SNR_FUN <- "+"
probs <- 0.5

configs <- expand.grid(SNR = SNR,
                       pi0 = pi0s,
                       SNR_FUN = SNR_FUN,
                       prob = probs)
seq_configs <- 1:nrow(configs)

# - - - - - - - - - - - - -
# experiments
# - - - - - - - - - - - - -
for (cc in seq_configs) {
    config <- configs[cc, ]
    
    pi0 <- config[["pi0"]]
    SNR <- config[["SNR"]]
    SNR_FUN <- config[["SNR_FUN"]]
    prob <- config[["prob"]]
    
    simname <- sprintf("%s_m=%s_pi0=%s_SNR=%s%s_prob=%s_B=%s_nb-exp=%s",
                       ds_name, m, pi0, SNR_FUN, SNR, prob,
                       B, nb_exp)
    print(simname)
    cat(cc, "/", nrow(configs), ":", simname, "\n")
    filename <- sprintf("%s.rds", simname)
    pathname <- file.path(path, filename)
    
    t0 <- Sys.time()
    # res <- lapply(1:nb_exp, FUN = function(i) {
    res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
        
        ## add some signal
        sig <- add_signal(X = X0, pi0 = pi0, SNR = SNR, SNR_FUN = SNR_FUN, prob = prob)
        
        ## define gene selections (for power estimation)
        tests <- rowTestFUN(sig$Y, sig$groups)
        p_values <- tests$p.value
        m <- length(p_values)
        rk <- rank(p_values)
        selections <- list(
            first_1 = which(rk %in% c(1:1)),
            first_10 = which(rk %in% c(1:10)),
            first_100 = which(rk %in% c(1:100)),
            first_1000 = which(rk %in% c(1:1000)),
            first_5000 = which(rk %in% c(1:5000)),
            first_10000 = which(rk %in% c(1:10000)),
            first_15000 = which(rk %in% c(1:15000)),
            BH_10 = which(p.adjust(p_values, method = "BH") < 0.10),
            BH_05 = which(p.adjust(p_values, method = "BH") < 0.05),
            p_05 = which(p_values < 0.05),
            p_01 = which(p_values < 0.01),
            H = 1:m)
        ## check JER control and estimate power
        res_i <- test_JER_control(
            Y = sig$Y, groups = sig$groups, truth = sig$truth, 
            rowTestFUN = rowWelchTests, B = B, 
            alpha = alphas, selections = selections,
            verbose = TRUE)
        list(level = tibble(exp = i, res_i$level),
             power = tibble(exp = i, res_i$power))
    })
    level <- Reduce(rbind, lapply(res, "[[", "level"))
    power <- Reduce(rbind, lapply(res, "[[", "power"))
    res <- list(level = level, power = power)
    saveRDS(res, file = pathname)
}
