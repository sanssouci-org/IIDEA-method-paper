library("sanssouci")
library("GSEABenchmarkeR")
library("future.apply")
library("cherry")
library("tidyr")
library("dplyr")
library("ggplot2")

# source("numerical-experiments/00_utils.R")

test_JER_control_DEG_intra <- function(Y, groups,  rowTestFUN = sanssouci::rowWelchTests, 
                                       B = 100, alpha = 0.05, 
                                       verbose = FALSE) {
  obj_i <- SansSouci(Y = Y, groups = groups)
  
  # Oracle predictions: depend on experiment (via gene order) but not on alpha!
  # obj_oracle <- obj_i
  # obj_oracle$input$truth <- truth
  # obj_oracle <- fit(obj_oracle, rowTestFUN = rowTestFUN, family = "Oracle", alpha = NA)
  # FP_oracle <- predict(obj_oracle, what = "FP", all = TRUE)$bound
  
  # TP_oracle <- sapply(selections, FUN = function(sel) {
  #   predict(obj_oracle, S = sel, what = "TP", all = FALSE)
  # })
  
  FDP0 <- TP0 <- NULL              ## Simes (no calibration)
  FDP <- TP <- NULL                ## Simes + single-step calibration
  FDP_sd <- TP_sd <- NULL          ## Simes + step-down calibration
  FDP_cherry <- TP_cherry <- NULL  ## cherry/ARI
  
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
    
    fdp <- predict(cal, what = "FDP", all = TRUE)
    # valid_bound <- all(FP >= FP_oracle)
    FDP[[aa]] <- tibble(alpha = alpha[aa], fdp)
    
    tp <- predict(cal,  what = "TP", all = TRUE)
    # pow <- format_power(TP / TP_oracle)
    TP[[aa]] <- tibble(alpha = alpha[aa], tp)
    
    ## Simes + calibration (step-down)
    cal_sd <- fit(cal, B = B, rowTestFUN = rowTestFUN, 
                  alpha = alpha[aa], family = "Simes")
    
    fdp_sd <- predict(cal_sd, what = "FDP", all = TRUE)
    # valid_bound <- all(FP >= FP_oracle)
    FDP_sd[[aa]] <- tibble(alpha = alpha[aa], fdp_sd)
    
    tp_sd <- predict(cal, what = "TP", all = TRUE)
    # pow_sd <- format_power(TP_sd / TP_oracle)
    TP_sd[[aa]] <- tibble(alpha = alpha[aa], tp_sd)
    
    ## Simes + no calibration
    thr <- t_linear(alpha[aa], 1:m, nHyp(cal)) ## hack to avoid recalc. pivotal stat
    cal0$output$thr <- thr
    
    fdp0 <- predict(cal0, what = "FDP", all = TRUE)
    # valid_bound <- all(FP0 >= FP_oracle)
    FDP0[[aa]] <- tibble(alpha = alpha[aa], fdp0)
    
    tp0 <- predict(cal0, what = "TP", all = TRUE)
    # pow0 <- format_power(TP0 / TP_oracle)
    TP0[[aa]] <- tibble(alpha = alpha[aa], tp0)
    
    ## cherry/ARI
    tp0_cherry <- curveSimes(hommel = hom, alpha = alpha[aa], plot = FALSE)
    fdp_cherry <- (1:m - tp0_cherry)/m
    # valid_bound <- all(FP_cherry >= FP_oracle)
    FDP_cherry[[aa]] <- tibble(alpha = alpha[aa], "bound" = fdp_cherry, "stat" = "FDP", x = 1:m, label = "cherry")
    
    # tp_cherry <- pickSimes(hommel = hom, alpha = alpha[aa], silent = TRUE)
    # pow_cherry <- format_power(TP_cherry / TP_oracle)
    TP_cherry[[aa]] <- tibble(alpha = alpha[aa], "bound"= tp0_cherry, "stat" = "TP", x = 1:m, label = "cherry")
    
  }
  FDP <- tibble(method = "Simes + single-step calibration", Reduce(rbind, FDP))
  FDP0 <- tibble(method = "Simes (parametric)", Reduce(rbind, FDP0))
  FDP_sd <- tibble(method = "Simes + step-down calibration", Reduce(rbind, FDP_sd))
  FDP_cherry <- tibble(method = "cherry/ARI", Reduce(rbind, FDP_cherry))
  
  TP <- tibble(method = "Simes + single-step calibration", Reduce(rbind, TP))
  TP0 <- tibble(method = "Simes (parametric)", Reduce(rbind, TP0))
  TP_sd <- tibble(method = "Simes + step-down calibration", Reduce(rbind, TP_sd))
  TP_cherry <- tibble(method = "cherry/ARI", Reduce(rbind, TP_cherry))
  
  list(FDP = rbind(FDP, FDP0, FDP_sd, FDP_cherry),
       TP = rbind(TP, TP0, TP_sd, TP_cherry))
}


# alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
alphas <- 0.1  # target JER level
B <- c(1e3)          # number of permutations for adaptive methods
nb_exp <- 1e3    # number of experiments

# future::availableCores() to know available 'workers'
plan(multisession, workers = 40)

technology <- c("microarray", "RNAseq")[2]

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
  
} else if (technology == "RNAseq") {
  rowTestFUN <- sanssouci::rowWilcoxonTests
  
  data("RNAseq_blca", package = "sanssouci.data")
  ds_name <- "BLCA"
  X <- RNAseq_blca
  groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
  rm(RNAseq_blca)
  
  # X0 <- X[, groups == 0]
  
  # filter out unexpressed genes
  BLCA0 <- X/colSums(X)*1e6
  ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
  if (length(ww) != 0){
    X <- X[-ww, ]
  }
  
}

# pi0 <- c(0.8)
# SNR = c(1.5, 2, 3)
# SNR_FUN = "*"
# prob <- 0.5

table(groups)
# groups0 <-  groups[groups == 0]

dim(X)
n <- dim(X)
m <- nrow(X)

Ns <- c(5,10,15,20,30, 40, 50, 60, 70, 80, 90, 100)

configs <- expand.grid(N = Ns)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_%s_sample_studies_FDP", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

######################## 02_run_differential_expression.R ###########

for (cc in seq_configs) {
  N <- configs[cc, "N"]
  # pi0 <- configs[cc, "pi0"]
  # SNR <- configs[cc, "SNR"]
  
  # X0_resize <- X0[,sample(1:length(groups0), N)]
  
  simname <- sprintf("%s_m=%s_N=%s_nb-exp=%s",
                     ds_name, m,  
                     N, nb_exp)
  print(simname)
  cat(cc, "/", nrow(configs), ":", simname, "\n")
  filename <- sprintf("%s.rds", simname)
  pathname <- file.path(path, filename)
  
  t0 <- Sys.time()
  # res <- lapply(1:nb_exp, FUN = function(i) {
  res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
    # res <- for(i in 1:nb_exp){
    
    groups0_resize = rep(0, N)
    # count = 0
    while(table(groups0_resize)[1] < 2 | table(groups0_resize)[1] > N-2){
      sampling <- sample(1:length(groups), N)
      X0_resize <- X[,sampling]
      if (length(which(rowSums(X0_resize) == 0))>0){
        X0_resize <- X0_resize[-which(rowSums(X0_resize) == 0),] ## be careful, here we change the #genes each experiments
      }
      groups0_resize <- groups[sampling]
      # count <- count +1
    }
    # count
    groups0_resize
    
    ## add some signal
    # sig <- add_signal(X = X0_resize, pi0 = pi0, SNR = SNR, SNR_FUN = SNR_FUN, prob = prob)
    # sig
    
    tests <- rowTestFUN(X0_resize, groups0_resize)
    p_values <- tests$p.value
    m <- length(p_values)
    rk <- rank(p_values)
    # selections <- list(
    #   # first_1 = which(rk %in% c(1:1)),
    #   first_10 = which(rk %in% c(1:10)),
    #   first_100 = which(rk %in% c(1:100)),
    #   first_1000 = which(rk %in% c(1:1000)),
    #   # first_5000 = which(rk %in% c(1:5000)),
    #   # first_10000 = which(rk %in% c(1:10000)),
    #   # first_15000 = which(rk %in% c(1:15000)),
    #   # BH_10 = which(p.adjust(p_values, method = "BH") < 0.10),
    #   BH_05 = which(p.adjust(p_values, method = "BH") < 0.05),
    #   p_05 = which(p_values < 0.05),
    #   # p_01 = which(p_values < 0.01),
    #   H = 1:m)
    ## check JER control and estimate power
    res_i <- test_JER_control_DEG_intra(
      Y = X0_resize, groups = groups0_resize,
      rowTestFUN = rowTestFUN, B = B, 
      alpha = alphas, 
      verbose = TRUE)
    list(FDP = tibble(exp = i, res_i$FDP),
         TP = tibble(exp = i, res_i$TP))
  })
  FDP <- Reduce(rbind, lapply(res, "[[", "FDP"))
  TP <- Reduce(rbind, lapply(res, "[[", "TP"))
  res <- list(FDP = FDP, TP = TP)
  saveRDS(res, file = pathname)
  print(simname)
}


######################"" 03_plot.R ###############

# technology <- c("microarray", "RNAseq")[2]
data_set <- switch(technology,
                   microarray = "GSE19188_m=21407",
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
  level <- tibble(level, pi0 = pi0, SNR = SNR, 
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
                  B = B, N=N, nb_exp = nb_exp)
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
  lev <- filter(level, pi0 == 0.8, 
                alpha <=0.75, N %in% c(5, 20, 50, 90), 
                SNR %in% c(1, 2))
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
       file = plotname, scale = 1, width = 16, height = 8)

## power
power <- Reduce(rbind, powList) 
dim(power)

power$method[power$method == "cherry/ARI"] <- "ARI"
power$method[power$method == "Simes (parametric)"] <- "Simes"
power$method[power$method == "Simes + step-down calibration"] <- "Adaptive Simes"
power$method[power$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

pow <- filter(power, 
              # estimate <= 1.0000001,
              # prob == 0.5,
              SNR ==2,
              selection %in% c("first_1000", "BH_05", "H"),
              alpha <= 0.9,
              N %in% c(5,20,50,90),
              pi0 == 0.8)
nb_exp_ <- unique(pow[["nb_exp"]])
pow

p <- ggplot(pow, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method, shape = method)) + 
  facet_grid(selection~N)+
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
