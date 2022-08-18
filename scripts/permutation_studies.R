library("sanssouci")
library("GSEABenchmarkeR")
library("future.apply")
library("cherry")
library("tidyr")
library("dplyr")
library("ggplot2")

source("00_utils.R")

alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
Bs <- c(100, 200, 500, 1000, 2000, 5000)          # number of permutations for adaptive methods
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
  BLCA0 <- X/colSums(X)*1e6
  ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
  if (length(ww) != 0){
    X <- X[-ww, ]
  }
  dim(X)
  
}

table(groups)

dim(X)
m <- nrow(X)

configs <- expand.grid(B = Bs)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_permutation_%s", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)

######################## 02_run_differential_expression.R ###########

for (cc in seq_configs) {
  B <- configs[cc, "B"]
  
  
  simname <- sprintf("%s_m=%s_B=%s_nb-exp=%s",
                     ds_name, m, 
                     B, nb_exp)
  print(simname)
  cat(cc, "/", nrow(configs), ":", simname, "\n")
  filename <- sprintf("%s.rds", simname)
  pathname <- file.path(path, filename)
  
  t0 <- Sys.time()
  # res <- lapply(1:nb_exp, FUN = function(i) {
  res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
    
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
      H = 1:m)
    ## check JER control and estimate power
    res_i <- test_JER_control_DEG_intra(
      Y = X, groups = groups, 
      rowTestFUN = rowTestFUN, B = B, 
      alpha = alphas, selections = selections,
      verbose = TRUE)
    list(lambda = tibble(exp = i, res_i$level),
         power = tibble(exp = i, res_i$power))
  })
  lambda <- Reduce(rbind, lapply(res, "[[", "lambda"))
  power <- Reduce(rbind, lapply(res, "[[", "power"))
  res <- list(lambda = lambda, power = power)
  saveRDS(res, file = pathname)
  print(simname)
}


######################"" 03_plot.R ###############

# technology <- c("microarray", "RNAseq")[2]
data_set <- switch(technology,
                   microarray = "GSE19188_m=21407",
                   RNAseq = "BLCA_m=12534")
path <- file.path("results", sprintf("diff-expr_permutation_%s", technology))
filenames <- list.files(path, pattern = data_set)
tail(filenames)
length(filenames)

pattern <- sprintf("%s_B=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)
filename <- filenames[1]

lambdaList <- list()
powList <- list()
for (filename in filenames) {
  print(filename)
  pathname <- file.path(path, filename)
  B <- as.numeric(gsub(pattern, "\\1", filename))
  nb_exp <- as.numeric(gsub(pattern, "\\2", filename))
  dat <- readRDS(pathname)
  
  lambda <- dat$lambda %>% group_by(alpha, method) 
  lambda <- tibble(lambda, B = B, nb_exp = nb_exp)
  lambdaList[[filename]] <- lambda
  
  power <- dat$power %>% group_by(alpha, method, selection) 
  power <- tibble(power, B = B, nb_exp = nb_exp)
  powList[[filename]] <- power
}

lambda <- Reduce(rbind, lambdaList) 
dim(lambda)

lambda$method[lambda$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

lev <- NULL

nb_exp_ <- unique(lev[["nb_exp"]])

sc_pct <- function(x, ...) scales::percent(x, accuracy = 1, ...)

# one big plot
p <- lambda %>%
  filter(alpha <= 0.5) %>%
  ggplot(aes(x = factor(B), y = lambda, group = B) ) +
  geom_boxplot(shape = "circle", fill = "#112446") +
  theme_minimal() +
  facet_wrap(vars(alpha)) + 
  ggtitle("Control of lambda") +
  xlab("Permutation B") + ylab("Lambda")
p
plotname <- sprintf("lambda-control_%s_%s_B=%s.pdf", technology, data_set, B)
ggsave(p, 
       file = plotname, scale = 1, width = 8, height = 8)

## power
power <- Reduce(rbind, powList) 
dim(power)

power$method[power$method == "cherry/ARI"] <- "ARI"
power$method[power$method == "Simes (parametric)"] <- "Simes"
power$method[power$method == "Simes + step-down calibration"] <- "Adaptive Simes"
power$method[power$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

p <- power %>%
  filter(alpha == 0.1, 
         selection %in% c("BH_05", "first_1000", "H")
  ) %>%
  mutate(B = factor(B), FDP = power) %>%
  ggplot() +
  aes(x = B, y = FDP, fill = B) +
  # geom_violin(adjust = 1L, scale = "area") +
  geom_violin() +
  theme_minimal() +
  facet_wrap(vars(selection), scale = "free") + 
  # ggtitle("Control of FDP by selection") + 
  scale_fill_brewer(palette = "Purples") + 
  #  scale_shape_discrete() +
  theme(legend.position="none")+
  stat_summary(fun = mean, geom="line", aes(group = 1)) + #
  # stat_summary(fun = median, geom="line", aes(group = 1)) + #
  stat_summary(fun = function(x){quantile(x, probs = 0.99)}, geom="line", aes(group = 1),linetype = "dashed") +
  stat_summary(fun = function(x){quantile(x, probs = 0.01)}, geom="line", aes(group = 1), linetype = "dashed") +
  xlab("Permutation (B)") + ylab("False Discoveries Proportion (for 1000 experiments)")
p
plotname <- sprintf("permutation_B-vs-bound_%s_%s.pdf", technology, data_set)
ggsave(p, file = plotname, scale = 1, width = 8, height = 6)

# power %>% filter(alpha ==0.05, selection == "first_1000") %>% mutate(B = factor(B), TP = power, FDP = 1-TP/1000)%>%
#   ggplot() +
#   aes(x = B, y = FDP) +
#   # geom_violin(adjust = 1L, scale = "area") +
#   geom_boxplot(shape = "circle") +
#   theme_minimal() +
#   facet_wrap(vars(selection), scales = "free") + 
#   ggtitle("Control of TP by selection") +
#   xlab("Permutation B") + ylab("True positive")

p <- power %>% 
  filter(selection %in% c("BH_05", "first_1000", "H")) %>%
  group_by(alpha, selection, B) %>% 
  summarise(mean_time = mean(time), 
            q_01 = quantile(power, probs = 0.01), 
            q_99 = quantile(power, probs = 0.99)) %>% 
  mutate(diff_q = q_99 - q_01) %>%
  filter(alpha == 0.1) %>%
ggplot(aes(x = mean_time, y = diff_q, color = selection, label = B)) + 
  geom_point() + 
  geom_line() +
  geom_text(hjust=-0.1, vjust=-0.1) +
  xlab("Mean time (across 1000 experiments)") + ylab("Interquantile range (Q0.99 - Q0.01) (across 1000 experiments)")
p
plotname <- sprintf("permutation_time-vs-IQR_%s_%s.pdf", technology, data_set)
ggsave(p, file = plotname, scale = 1, width = 8, height = 6)
p1 <- power %>% select(B, time) %>% mutate(B = factor(B)) %>% group_by(B) %>% ggplot() + geom_violin(aes(y = time, x = B)) +
  xlab("Permutation (B)") + ylab("Time (across 1000 experiments)")
p1
ggsave(p1, file = "boxplot_time_permutation.pdf", scale = 1, width = 8, height = 6)

power %>% 
  filter(selection %in% c("BH_05", "first_1000", "H")) %>%
  group_by(alpha, selection, B) %>% 
  summarise(mean_time = mean(time)) %>% filter(B == 1000) %>% ungroup() %>% select(mean_time) %>% unique()

power %>% filter(B == 1000, alpha == 0.1, selection == "BH_05") %>%
  group_by(B) %>% 
  summarise(q0.995 = quantile(power, probs = 0.995), 
            q0.005 = quantile(power, probs = 0.005), 
            q0.99 = quantile(power, probs = 0.99),
            q0.01 = quantile(power, probs = 0.01))







p <- power %>%
  filter(alpha == 0.1, 
         selection %in% c("BH_05", "first_1000", "H")
  ) %>%
  # mutate(B = factor(B)) %>% 
  mutate(FDP = power)%>% 
  left_join(power %>% group_by(B) %>% summarise(t = mean(time)), by="B")%>%
  ggplot() +
  aes(x = (as.numeric(t)), y = FDP, fill = factor(B), group=t) +
  # geom_violin(adjust = 1L, scale = "area") +
  geom_violin(scale = "width") +
  theme_minimal() +
  facet_wrap(vars(selection), scale = "free") + 
  # ggtitle("Control of FDP by selection") + 
  scale_fill_brewer(palette = "Purples") + 
  #  scale_shape_discrete() +
  theme(legend.position="bottom")+
  stat_summary(fun = mean, geom="line", aes(group = 1)) + #
  scale_x_continuous(trans = "log10", breaks = c(10, 15, 20, 30, 50))  +
  # stat_summary(fun = median, geom="line", aes(group = 1)) + #
  stat_summary(fun = function(x){quantile(x, probs = 0.99)}, geom="line", aes(group = 1),linetype = "dashed") +
  stat_summary(fun = function(x){quantile(x, probs = 0.01)}, geom="line", aes(group = 1), linetype = "dashed") +
  xlab("Mean time in log10 scale (sec)") + ylab("False Discoveries Proportion (across 1000 experiments)") + labs(fill = "Permutation (B)")
p
ggsave(p, file = "permutation_violin_time_FDP_BLCA_m=12534.pdf", scale = 1, width = 8, height = 6)
