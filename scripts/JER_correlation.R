library("sanssouci")
# library("GSEABenchmarkeR")
library("future.apply")
library("cherry")
library("tidyr")
library("dplyr")
library("ggplot2")
library(matrixStats)

format_power <- function(pow) {
  as.data.frame(pow) %>% 
    tibble::rownames_to_column(var = "selection") %>%
    tidyr::pivot_longer(cols = !`selection`, values_to = "power") %>%
    dplyr::select(!name) %>% 
    filter(is.finite(power))  ## only keep those entries with non 0 oracle TP (ie condition on |S \cap H1|>0)
}

test_JER_control_DEG_intra <- function(Y, groups, truth, rowTestFUN = sanssouci::rowWelchTests, 
                                       B = 100, alpha = 0.05, 
                                       selections = NULL, verbose = FALSE) {
  # obj_i <- SansSouci(Y = X, groups = groups)
  
  # Oracle predictions: depend on experiment (via gene order) but not on alpha!
  # obj_oracle <- obj_i
  # obj_oracle$input$truth <- truth
  
  # obj_oracle <- fit(obj_oracle, rowTestFUN = rowTestFUN, family = "Oracle", alpha = NA)
  p.value <- rowTestFUN(Y, groups)$p.value
  # FP_oracle <- predict(obj_oracle, what = "FP", all = TRUE)$bound
  FP_oracle <- sanssouci:::posthoc_bound(p.values = p.value, thr = truth, lab = "Simes", all = TRUE, what = "FP")$bound
  
  TP_oracle <- sapply(selections, FUN = function(sel) {
    # predict(obj_oracle, S = sel, what = "TP", all = FALSE)
    sanssouci:::posthoc_bound(p.values = p.value, S = sel, thr = truth, lab = "Simes", all = FALSE, what = "TP")$bound
  })
  
  level0 <- power0 <- NULL              ## Simes (no calibration)
  level <- power <- NULL                ## Simes + single-step calibration
  level_sd <- power_sd <- NULL          ## Simes + step-down calibration
  level_cherry <- power_cherry <- NULL  ## cherry/ARI
  
  # cal0 <- fit(obj_i, B = 0, rowTestFUN = rowTestFUN, 
  #             alpha = alpha[1],
  #             family = "Simes")
  perm <- get_perm(Y, groups, B, rowTestFUN = cor_test2, alternative = "two.sided")
  perm_p_value <- perm$p.value
  
  calibration0 <- calibrate(p0 = perm_p_value, m = nrow(X), alpha = alpha[1])
  
  # cal <- fit(obj_i, B = B, rowTestFUN = rowTestFUN, 
  #            alpha = alpha[1],
  #            family = "Simes", max_steps_down = 0)
  calibration <- calibrate(p0 = perm_p_value, m = nrow(X), alpha = alpha[1], max_steps_down = 0)
  
  piv_stat <- calibration$piv_stat  # does not depend on alpha[1] because single-step
  m <- dim(Y)[1]
  
  p_values <- p.value
  hom <- hommelFast(p_values)
  
  for (aa in seq(along = alpha)) {
    ## Simes + calibration (single-step)
    lambda <- stats::quantile(piv_stat, alpha[aa], type = 1)
    thr <- t_linear(lambda, 1:m, dim(Y)[1]) ## hack to avoid recalc. pivotal stat
    # cal$output$thr <- thr
    
    # FP <- predict(cal, what = "FP", all = TRUE)$bound
    FP <- sanssouci:::posthoc_bound(p.values = p.value, thr = thr, lab = "Simes", all = TRUE, what = "FP")$bound
    valid_bound <- all(FP >= FP_oracle)
    level[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)
    
    TP <- sapply(selections, FUN = function(sel) {
      # predict(cal, S = sel, what = "TP", all = FALSE)
      sanssouci:::posthoc_bound(p.values = p.value, S = sel, thr = thr, lab = "Simes", all = FALSE, what = "TP")$bound
    })
    pow <- format_power(TP / TP_oracle)
    power[[aa]] <- tibble(alpha = alpha[aa], pow)
    
    ## Simes + calibration (step-down)
    # cal_sd <- fit(cal, B = B, rowTestFUN = rowTestFUN, 
    #               alpha = alpha[aa], family = "Simes")
    calibration_sd <- calibrate(p0 = perm_p_value, m = nrow(X), alpha = alpha[aa])
    
    # FP <- predict(cal_sd, what = "FP", all = TRUE)$bound
    FP <- sanssouci:::posthoc_bound(p.values = p.value, thr = calibration_sd$thr, lab = "Simes", all = TRUE, what = "FP")$bound
    valid_bound <- all(FP >= FP_oracle)
    level_sd[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)
    
    TP_sd <- sapply(selections, FUN = function(sel) {
      # predict(cal, S = sel, what = "TP", all = FALSE)
      sanssouci:::posthoc_bound(p.values = p.value, S = sel, thr = calibration_sd$thr, lab = "Simes", all = FALSE, what = "TP")$bound
    })
    pow_sd <- format_power(TP_sd / TP_oracle)
    power_sd[[aa]] <- tibble(alpha = alpha[aa], pow_sd)
    
    ## Simes + no calibration
    thr <- t_linear(alpha[aa], 1:m, dim(Y)[1]) ## hack to avoid recalc. pivotal stat
    
    # FP0 <- predict(cal0, what = "FP", all = TRUE)$bound
    FP0 <- sanssouci:::posthoc_bound(p.values = p.value, thr = thr, lab = "Simes", all = TRUE, what = "FP")$bound
    valid_bound <- all(FP0 >= FP_oracle)
    level0[[aa]] <- tibble(alpha = alpha[aa], "valid bound" = valid_bound)
    
    TP0 <- sapply(selections, FUN = function(sel) {
      # predict(cal0, S = sel, what = "TP", all = FALSE)
      sanssouci:::posthoc_bound(p.values = p.value, S = sel, thr = thr, lab = "Simes", all = FALSE, what = "TP")$bound
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
  
  list(level = rbind(level, level0, level_sd, level_cherry),
       power = rbind(power, power0, power_sd, power_cherry))
}



alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
B <- c(1e3)          # number of permutations for adaptive methods
nb_exp <- 1e3

# future::availableCores() to know available 'workers'
plan(multisession, workers = 40)

technology <- c("microarray", "RNAseq", "simulation")[3]
if (technology == "microarray") {
  
  data("expr_ALL", package = "sanssouci.data")
  X <- expr_ALL
  groups <- colnames(X)
  groups[which(groups == "BCR/ABL")] <- 1
  groups[which(groups == "NEG")] <- 0
  # geo2kegg <-  R.cache::memoizedCall(loadEData,"geo2kegg")
  ds_name <- "leukemia"
  # rawData <- R.cache::memoizedCall(maPreproc, geo2kegg[ds_name])[[1]]
  # X <- SummarizedExperiment::assays(rawData)$exprs
  # cats <- SummarizedExperiment::colData(rawData)
  # ww <- match(cats$Sample, base::colnames(X))
  # groups <- cats$GROUP[ww]
  
} else if (technology == "RNAseq") {
  
  data("RNAseq_blca", package = "sanssouci.data")
  ds_name <- "BLCA"
  X <- RNAseq_blca
  groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
  rm(RNAseq_blca)
  
  X0 <- X[, groups == 0]
  
  # filter out unexpressed genes
  BLCA0 <- X0/colSums(X)*1e6
  ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
  if (length(ww) != 0){
    X0 <- X0[-ww, ]
  }
  
} else if (technology == "simulation"){
  n <- 60
  p <- 10000
  rho <- c(0, 0.1, 0.4, 0.6, 0.8, 1)
  K <- 99
  
  ds_name = "simulation"
  
  
  mu <- rep(0, p)
  # Sigma <- toeplitz(rho ^ (0:(p-1)))[1:5, 1:5]
  
  truthFp <- c(1, sort(sample(1:p, K, replace=F)), p)
  # truthS <- sort(sample(1:K, p, replace = T))
  
  
  truthF <- rep(0, p)
  for(i in 1:(length(truthFp)-1)){
    indF <- truthFp[i]:truthFp[i+1]
    truthF[indF] <- i
    # Sigma[indF, indF] <- Sigma[indF,indF] + rho
  }
  # for (i in unique(truthF)){
  #   # indS <- which(truthS == i)
  #   indF <- which(truthF == i)
  #   Sigma[indF, indF] <- Sigma[indF,indF] + rho
  # }
  # diag(Sigma) <- 1
  # # image(Sigma)
  # X <- MASS::mvrnorm(n= n, mu = mu, Sigma = Sigma)
  # X <- t(X)
  
}


#### Y 
# Y <- rnorm(dim(X)[2], 0, 1)
# Y <- rbinom(dim(X)[2], prob = 0.1, size = 10)
# Y = X[1,]*2+X[2,]*3
# sel = 1:50
# Y <- colSums(X[sel,])


cor_test2 <- function(X,categ,alternative = "two.sided"){
  r <- matrixTests::row_cor_pearson(X, categ, alternative = alternative)
  return(list(p.value = r$pvalue, statistic = r$statistic, estimate =r$cor))
}

rowTestFUN <- cor_test2

# dim(X)
m <- p

snr <- c(2)
effect <- 1



# nb_correlation = c(0.01,0.05, 0.1, 0.15,0.2,0.25,0.5)
nb_correlation = c(1, 2, 5, 10, 20, 50)#1:length(unique(truthF))

configs <- expand.grid(nb_correlation = nb_correlation)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_%s_correlation", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)
for (r in rho){
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (i in unique(truthF)){
    # indS <- which(truthS == i)
    indF <- which(truthF == i)
    Sigma[indF, indF] <- Sigma[indF,indF] + r
  }
  diag(Sigma) <- 1
  # image(Sigma)
  X <- MASS::mvrnorm(n= n, mu = mu, Sigma = Sigma)
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
      # res <- for(i in 1:nb_exp){
      
      ## add some signal
      # truth <- rep(0, dim(X)[1])
      # truth[sample(1:length(truth), nb_cor)] <- 1
      # Y <- truth %*% X +  rnorm(dim(X)[2], 0, 1)
      
      select <- sample(unique(truthF), nb_cor)
      non_zero = which(truthF %in% select)
      truth = rep(0, p)
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
      res_i <- test_JER_control_DEG_intra(
        Y = X, groups = Y, truth = truth, 
        rowTestFUN = rowTestFUN, B = B, 
        alpha = alphas, selections = selections,
        verbose = TRUE)
      len <- sapply(selections, FUN = length) %>% data.frame() %>% tibble::rownames_to_column("selection") %>% rename('length' = '.')
      list(level = tibble(exp = i, res_i$level),
           power = tibble(exp = i, res_i$power) %>% left_join(., len, by = "selection"))
    })
    level <- Reduce(rbind, lapply(res, "[[", "level"))
    power <- Reduce(rbind, lapply(res, "[[", "power"))
    res <- list(level = level, power = power)
    saveRDS(res, file = pathname)
    print(simname)
  }
}


######################"" 03_plot.R ###############

# technology <- c("microarray", "RNAseq", "simulation")[3]
data_set <- switch(technology,
                   microarray = "leukemia_m=9038",
                   RNAseq = "BLCA_m=12418", 
                   simulation = "simulation_m=10000")
path <- file.path("results", sprintf("diff-expr_%s_correlation", technology))
filenames <- list.files(path, pattern = data_set)
tail(filenames)
length(filenames)

pattern <- sprintf("%s_nb_cor=(.*)_rho=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)
filename <- filenames[1]

levList <- list()
powList <- list()
for (filename in filenames) {
  print(filename)
  pathname <- file.path(path, filename)
  # nb_cor <- as.numeric(gsub(pattern, "\\1", filename))
  nb_cor <- as.numeric(gsub(pattern, "\\1", filename))
  rho <- round(as.numeric(gsub(pattern, "\\2", filename)), 1)
  nb_exp <- as.numeric(gsub(pattern, "\\3", filename))
  dat <- readRDS(pathname)
  
  level <- dat$level %>% group_by(alpha, method) %>%
    summarise(estimate = 1 - mean(`valid bound`), 
              std_err = sqrt((1 - estimate)*estimate/n()),
              min = estimate - 2*std_err,
              max = estimate + 2*std_err,
              .groups = "drop")
  level <- tibble(level, nb_cor = nb_cor,
                  B = B, nb_exp = nb_exp, rho = rho)
  levList[[filename]] <- level
  
  power <- dat$power %>% group_by(alpha, method, selection) %>%
    summarise(estimate = mean(power, na.rm = TRUE),
              n = sum(!is.na(power)), # to check if large enough support for estimation
              std_err = sd(power, na.rm = TRUE)/sqrt(n),
              min = estimate - 2*std_err,
              max = estimate + 2*std_err,
              length = mean(length),
              .groups = "drop")
  power <- tibble(power, nb_cor = nb_cor,
                  B = B, nb_exp = nb_exp, rho = rho)
  powList[[filename]] <- power
}

level <- Reduce(rbind, levList) 
dim(level)

level$method[level$method == "cherry/ARI"] <- "ARI"
level$method[level$method == "Simes (parametric)"] <- "Simes"
level$method[level$method == "Simes + step-down calibration"] <- "Adaptive Simes"
level$method[level$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

lev <- level
if (technology == "microarray") {
  lev <- filter(level, alpha <= 0.5, pi0 %in% c(0.5, 0.8, 0.9, 1), prob == 0.5)
} else {
  lev <- filter(level, alpha <=0.7, 
                nb_cor %in% c(1,5,20),
                rho %in% c(0.2, 0.4, 0.6),
                # method %in% c("ARI", "Adaptive Simes")
                )
}

nb_exp_ <- unique(lev[["nb_exp"]])

sc_pct <- function(x, ...) scales::percent(x, accuracy = 1, ...)

# one big plot
p <- ggplot(lev, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method,
                     fill = method,
                     shape = method)) + 
  # facet_grid(SNR~N, labeller = label_bquote(cols = N: .(N),
  #                                           rows = SNR: .(SNR),))  +
  facet_grid(rho~nb_cor, labeller = label_bquote(cols = pi[0]: .(1-nb_cor/100),
                                                 rows = rho: .(rho))) +
  geom_point() + 
  geom_line() + 
  geom_abline(slope = 1, intercept = 0) + 
  ylab(paste("Empirical risk achieved (estimated from", nb_exp_, "experiments)")) +
  xlab(expression(paste("Target risk (", alpha, ")"))) +
  labs(color = "Method", fill = "Method", shape = "Method")+
  scale_x_continuous(labels = sc_pct) + 
  scale_y_continuous(labels = sc_pct) + 
  theme(legend.position = "bottom")
p + geom_ribbon(alpha = 0.3, linetype = 1)
plotname <- sprintf("JER-control_%s_%s_correlation.pdf", technology, data_set)
ggsave(p + geom_ribbon(alpha = 0.3, linetype = 1), 
       file = plotname, scale = 1, width = 8, height = 8)

## power
power <- Reduce(rbind, powList) 
dim(power)

power$method[power$method == "cherry/ARI"] <- "ARI"
power$method[power$method == "Simes (parametric)"] <- "Simes"
power$method[power$method == "Simes + step-down calibration"] <- "Adaptive Simes"
power$method[power$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

pow <- power
pow <- filter(power,
              # estimate <= 1.0000001,
              # prob == 0.5,
              alpha <=0.5, 
              nb_cor %in% c(1,5,20),
              rho %in% c(0.4),
              selection %in% c("BH_05", "first_1000", "H"),
)
nb_exp_ <- unique(pow[["nb_exp"]])
pow

p <- ggplot(pow, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method, shape = method)) + 
  facet_grid(selection~nb_cor, labeller = label_bquote(cols = pi[0]: .(1-nb_cor/100),
                                                       ))+
  geom_point() + 
  geom_line() + 
  ylim(c(0, 1.0001)) +
  ylab(paste("Power (estimated from", nb_exp_, "experiments)")) +
  xlab(expression(paste("Target risk (", alpha, ")"))) +
  labs(color = "Method", fill = "Method", shape = "Method")+
  scale_x_continuous(labels = sc_pct) + 
  scale_y_continuous(labels = sc_pct) + 
  theme(legend.position="bottom")
p
plotname <- sprintf("power_%s_%s_correlation.pdf", technology, data_set)
ggsave(p, file = plotname, scale = 1, width = 8, height = 6)
