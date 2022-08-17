library("sanssouci")
library("GSEABenchmarkeR")
library("future.apply")
library("cherry")
library("tidyr")
library("dplyr")
library("ggplot2")

source("numerical-experiments/00_utils.R")

test_JER_control_DEG_intra <- function(Y, groups, truth, rowTestFUN = sanssouci::rowWelchTests, 
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
  
  level0 <- power0 <- NULL              ## Simes (no calibration)
  level <- power <- NULL                ## Simes + single-step calibration
  level_sd <- power_sd <- NULL          ## Simes + step-down calibration
  level_cherry <- power_cherry <- NULL  ## cherry/ARI
  
  cal0 <- fit(obj_i, B = 0, rowTestFUN = rowTestFUN, 
              alpha = alpha[1],
              family = "Simes")
  
  cal <- fit(obj_i, B = B, rowTestFUN = rowTestFUN, 
             alpha = alpha[1],
             family = "Simes", max_steps_down = 0)
  piv_stat <- cal$output$piv_stat  # does not depend on alpha[1] because single-step
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
    cal_sd <- fit(cal, B = B, rowTestFUN = rowTestFUN, 
                  alpha = alpha[aa], family = "Simes")
    
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
  
  list(level = rbind(level, level0, level_sd, level_cherry),
       power = rbind(power, power0, power_sd, power_cherry))
}

# rowWilcoxonTests <- function (mat, categ, alternative = c("two.sided", "less", "greater"), 
#                               correct = TRUE) 
# {
#   alternative <- match.arg(alternative)
#   stopifnot(all(categ %in% c(0, 1)))
#   categ <- as.matrix(categ)
#   levels(categ) <- NULL
#   # apply(categ, 2, categCheck, n = ncol(mat))
#   B <- ncol(categ)
#   m <- nrow(mat)
#   n_obs <- nrow(categ)
#   
#   rks <- rowRanks(mat, ties.method = "average")
#   nx <- colSums(categ)
#   ny <- n_obs - nx
#   
#   min_stat <- nx * (nx + 1)/2 
#   
#   stats <- rks %*% categ
#   # stats <- stats - min_stat 
#   stats <- sweep(stats, MARGIN = 2, STATS = min_stat, FUN = "-")
#   
#   rks2 <- rks 
#   mode(rks2) <- "integer"
#   n_ties <- rowTabulates(rks2) 
#   ties <- rowSums(n_ties^3 - n_ties) 
#   # sigma <- sqrt((nx * ny/12) * ((ny + nx + 1) - ties/((ny + nx) * (ny + nx - 1))))
#   # sigma <- t(sqrt((nx * ny/12) * ((ny + nx + 1) - matrix(rep(ties, B), nrow = B, byrow = T)/((ny + nx) * (ny + nx - 1)))))
#   quotient <- sweep(matrix(rep(ties, B), ncol = B), MARGIN = 2, STATS = (ny + nx) * (ny + nx - 1), FUN = "/")
#   difference <- sweep(-quotient, MARGIN = 2, STATS = (ny + nx + 1), FUN = "+")
#   sigma <- sqrt(sweep(difference, MARGIN = 2, STATS = (nx * ny/12), FUN = "*"))
#   
#   # z <- stats - ny * nx/2
#   z <- sweep(stats, MARGIN = 2, STATS = ny * nx/2, FUN = "-")
#   CORRECTION <- 0
#   if (correct) {
#     CORRECTION <- switch(alternative, two.sided = sign(z) * 
#                            0.5, greater = 0.5, less = -0.5)
#   }
#   
#   z <- (z - CORRECTION)/sigma
#   p <- switch(alternative, 
#               less = pnorm(z), 
#               greater = pnorm(z, lower.tail = FALSE), 
#               two.sided = 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE)))
#   p[which(p<0)] <- 0
#   p[which(p>1)] <- 1
#   
#   est <- matrix(NA_real_, nrow = m, ncol = B)
#   if (dim(categ)[2] == 1){
#     wx <- which(categ == 1)
#     est <- rowMedians(mat[, wx]) - rowMedians(mat[, -wx])
#     stats <- as.vector(stats)
#     p <- as.vector(p)
#   } 
#   list(p.value = p, statistic = stats, estimate = est)
# }


alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
B <- c(1e3)          # number of permutations for adaptive methods
nb_exp <- 1e3    # number of experiments

# future::availableCores() to know available 'workers'
plan(multisession, workers = 40)

technology <- c("microarray", "RNAseq")[1]

# create data set and experiment parameters
if (technology == "microarray") {
  rowTestFUN <- sanssouci::rowWelchTests
  
  geo2kegg <-  R.cache::memoizedCall(loadEData,"geo2kegg")
  ds_name <- "GSE19188"
  rawData <- R.cache::memoizedCall(maPreproc, geo2kegg[ds_name])[[1]]
  X <- SummarizedExperiment::assays(rawData)$exprs
  cats <- SummarizedExperiment::colData(rawData)
  ww <- match(cats$Sample, base::colnames(X))
  groups <- cats$GROUP[ww]
  
  X0 <- X[, groups == 1]
  groups0 <-  groups[groups == 1]
  
} else if (technology == "RNAseq") {
  rowTestFUN <- sanssouci::rowWilcoxonTests
  
  data("RNAseq_blca", package = "sanssouci.data")
  ds_name <- "BLCA"
  X <- RNAseq_blca
  groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
  rm(RNAseq_blca)
  
  X0 <- X[, groups == 0]
  groups0 <-  groups[groups == 0]
  
  # filter out unexpressed genes
  BLCA0 <- X0/colSums(X)*1e6
  ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
  if (length(ww) != 0){
    X0 <- X0[-ww, ]
  }
  
}

pi0 <- c(0.8)
SNR = c(1, 1.5, 2, 3)
SNR_FUN = "*"
prob <- 0.5

table(groups)


dim(X0)
m <- nrow(X0)

Ns <- c(8,9)

configs <- expand.grid(N = Ns, pi0 = pi0, SNR = SNR)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_%s_sample_studies", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

######################## 02_run_differential_expression.R ###########

for (cc in seq_configs) {
  N <- configs[cc, "N"]
  pi0 <- configs[cc, "pi0"]
  SNR <- configs[cc, "SNR"]
  
  # X0_resize <- X0[,sample(1:length(groups0), N)]
  
  simname <- sprintf("%s_m=%s_SNR=%s_pi0=%s_N=%s_nb-exp=%s",
                     ds_name, m, SNR, pi0, 
                     N, nb_exp)
  print(simname)
  cat(cc, "/", nrow(configs), ":", simname, "\n")
  filename <- sprintf("%s.rds", simname)
  pathname <- file.path(path, filename)
  
  t0 <- Sys.time()
  # res <- lapply(1:nb_exp, FUN = function(i) {
  res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
  # res <- for(i in 1:nb_exp){
    s <- sample(1:length(groups0), N)
    X0_resize <- X0[,s]
    if (length(which(rowSums(X0_resize) == 0))>0){
      X0_resize <- X0_resize[-which(rowSums(X0_resize) == 0),] ## be careful, here we change the #genes each experiments
    }
    
    
    ## add some signal
    sig <- add_signal(X = X0_resize, pi0 = pi0, SNR = SNR, SNR_FUN = SNR_FUN, prob = prob)
    if(N <= 10){
      null_groups2 <- t(gtools::permutations(n = 2, r = length(sig$groups), v = c(0,1), repeats.allowed = TRUE))
      nonsame <- which(rowSums((sig$Y%*%(-(null_groups2-1)) == sig$Y%*%null_groups2)*1) ==0)
      sig$Y <- sig$Y[nonsame,]
      sig$truth <- sig$truth[nonsame]
    }
    noconstant <- which(rowVars(sig$Y)!=0)
    sig$Y <- sig$Y[noconstant,]
    sig$truth <- sig$truth[noconstant]
    
    
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
    res_i <- test_JER_control_DEG_intra(
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
  print(simname)
}


######################"" 03_plot.R ###############

# technology <- c("microarray", "RNAseq")[1]
data_set <- switch(technology,
                   microarray = "GSE19188_m=21405",
                   RNAseq = "BLCA_m=12418")
path <- file.path("results", sprintf("diff-expr_%s_sample_studies", technology))
filenames <- list.files(path, pattern = data_set)
tail(filenames)
length(filenames)

pattern <- sprintf("%s_SNR=(.*)_pi0=(.*)_N=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)
filename <- filenames[1]

levList <- list()
powList <- list()
for (filename in filenames) {
  print(filename)
  pathname <- file.path(path, filename)
  SNR <- as.numeric(gsub(pattern, "\\1", filename))
  pi0 <- as.numeric(gsub(pattern, "\\2", filename))
  N <- as.numeric(gsub(pattern, "\\3", filename))
  nb_exp <- as.numeric(gsub(pattern, "\\4", filename))
  dat <- readRDS(pathname)
  
  level <- dat$level %>% group_by(alpha, method) %>%
    summarise(estimate = 1 - mean(`valid bound`), 
              std_err = sqrt((1 - estimate)*estimate/n()),
              min = estimate - 2*std_err,
              max = estimate + 2*std_err,
              .groups = "drop")
  level <- tibble(level, pi0 = pi0, SNR = SNR, prob = prob,
                  B = B, N=N, nb_exp = nb_exp)
  levList[[filename]] <- level
  
  power <- dat$power %>% group_by(alpha, method, selection) %>%
    summarise(estimate = mean(power, na.rm = TRUE),
              n = sum(!is.na(power)), # to check if large enough support for estimation
              std_err = sd(power, na.rm = TRUE)/sqrt(n),
              min = estimate - 2*std_err,
              max = estimate + 2*std_err,
              .groups = "drop")
  power <- tibble(power, pi0 = pi0, SNR = SNR, 
                  prob = prob, B = B, N=N, nb_exp = nb_exp)
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
  lev <- filter(level, alpha <= 0.7, pi0 %in% c(0.5, 0.8, 0.9, 1), prob == 0.5, N %in% c(10, 50, 90))
} else {
  lev <- filter(level,pi0 == 0.8, alpha = 0.7 )
}

nb_exp_ <- unique(lev[["nb_exp"]])

sc_pct <- function(x, ...) scales::percent(x, accuracy = 1, ...)

# one big plot
p <- ggplot(lev, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method,
                     fill = method,
                     shape = method)) + 
  facet_grid(SNR~N, labeller = label_bquote(cols = N: .(N),
                                              rows = SNR: .(SNR),))  +
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
plotname <- sprintf("JER-control_sample_studies_%s_%s.pdf", technology, data_set)
ggsave(p + geom_ribbon(alpha = 0.3, linetype = 1), 
       file = plotname, scale = 1, width = 8, height = 8)

## power
power <- Reduce(rbind, powList) 
dim(power)

power$method[power$method == "cherry/ARI"] <- "ARI"
power$method[power$method == "Simes (parametric)"] <- "Simes"
power$method[power$method == "Simes + step-down calibration"] <- "Adaptive Simes"
power$method[power$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

pow <- filter(power, 
              estimate <= 1.0000001,
              # prob == 0.5,
              # SNR %in% c("*2"),
              # selection %in% c("first_100", "BH_05", "p_05", "H"), 
              selection == "H",
              alpha <= 0.7,
              pi0 == 0.8)
nb_exp_ <- unique(pow[["nb_exp"]])
pow

p <- ggplot(pow, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method, shape = method)) + 
  facet_grid(SNR~N, labeller = label_bquote(cols = N: .(N),
                                            rows = SNR: .(SNR),))+
  geom_point() + 
  geom_line() + 
  # ylim(c(0, 1.0001)) +
  ylab(paste("Power (estimated from", nb_exp_, "experiments)")) +
  xlab(expression(paste("Target risk (", alpha, ")"))) +
  labs(color = "Method", fill = "Method", shape = "Method")+
  scale_x_continuous(labels = sc_pct) + 
  scale_y_continuous(labels = sc_pct) + 
  theme(legend.position="bottom")
p
plotname <- sprintf("power_sample_studies_%s_%s.pdf", technology, data_set)
ggsave(p, file = plotname, scale = 1, width = 8, height = 6)
