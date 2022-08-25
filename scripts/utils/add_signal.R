add_signal <- function(
    X, pi0 = 0.7, SNR = 1, SNR_FUN = "+", prob = 0.5, verbose = TRUE) {
    n_genes <- nrow(X)
    n_obs <- ncol(X)
    
    # groups of observations
    groups <- rbinom(n_obs, 1, prob)
    while (sum(groups) <= 1 || sum(groups) >= n_obs-1) {
        print(groups)
        groups <- rbinom(n_obs, 1, prob)
    }    
    
    # differentially expressed genes
    n_DE <- round(n_genes*(1 - pi0))
    DE <- sample(n_genes, n_DE)
    truth <- rep(0, n_genes)
    truth[DE] <- 1
    stopifnot(all(truth[DE] == 1))
    stopifnot(all(truth[-DE] == 0))
    
    # adding signal for DE genes
    Y <- X
    if (SNR_FUN == "+") {
        delta <- runif(n_DE, min = -1, max = 1)*SNR/sqrt(n_obs)
        Y[DE, groups == 1] <- X[DE, groups == 1] + delta
    } else if (SNR_FUN == "*") {
        z <- SNR
        delta <- runif(n_DE, min = 1, max = z)
        b <- rbinom(n_DE, 1, 0.5)
        Y[DE[b == 1], groups == 1] <- X[DE[b == 1], groups == 1] * delta[b == 1]
        Y[DE[b == 0], groups == 0] <- X[DE[b == 0], groups == 0] * delta[b == 0]
    } else {
        stop("Argument 'SNR_FUN' should be '+' or '*'")
    }
    
    list(Y = Y, groups = groups, truth = truth)
}
