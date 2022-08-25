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
# parallelization setup
# - - - - - - - - - - - - -
plan(multisession, workers = 2) # by default parallelize on 2 nodes
# future::availableCores() to know available 'workers'

# - - - - - - - - - - - - -
# data set
# - - - - - - - - - - - - -
technology <- "RNAseq"
rowTestFUN <- sanssouci::rowWilcoxonTests
ds_name <- "BLCA"

source("scripts/utils/load_RNAseq_data.R") # loads 'X0' and 'groups'
m <- nrow(X0)
data_set <- sprintf("%s_m=%s", ds_name, m)
print(data_set)

path <- sprintf("results/diff-expr_%s_sample-size", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
B <- 10          # number of permutations for adaptive methods
nb_exp <- 10     # number of experiments

# Ns <- c(10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)
Ns <- c(10, 50, 90)
pi0 <- c(0.8)
# SNR <- c(1, 1.5, 2, 3)
SNR <-  c(1, 2)
SNR_FUN = "*"
prob <- 0.5

# - - - - - - - - - - - - -
# experiments
# - - - - - - - - - - - - -
configs <- expand.grid(N = Ns, pi0 = pi0, SNR = SNR)
seq_configs <- 1:nrow(configs)

for (cc in seq_configs) {
    N <- configs[cc, "N"]
    pi0 <- configs[cc, "pi0"]
    SNR <- configs[cc, "SNR"]
    n0 <- ncol(X0)
    
    simname <- sprintf("%s_SNR=%s_pi0=%s_N=%s_nb-exp=%s",
                       data_set, SNR, pi0, N, nb_exp)
    print(simname)
    cat(cc, "/", nrow(configs), ":", simname, "\n")
    filename <- sprintf("%s.rds", simname)
    pathname <- file.path(path, filename)
    
    t0 <- Sys.time()
    res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
        X0_resize <- X0[, sample(1:n0, N)]
        null_counts <- which(rowSums(X0_resize) == 0)
        if (length(null_counts)>0){
            ## caution: here we change the #genes in each experiment
            X0_resize <- X0_resize[-null_counts, ] 
        }
        
        ## add some signal
        sig <- add_signal(X = X0_resize, pi0 = pi0, 
                          SNR = SNR, SNR_FUN = SNR_FUN, 
                          prob = prob)
        str(sig)
        
        tests <- rowTestFUN(sig$Y, sig$groups)
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
        res_i <- test_JER_control(
            Y = sig$Y, groups = sig$groups, truth = sig$truth, 
            rowTestFUN = rowTestFUN, B = B, 
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