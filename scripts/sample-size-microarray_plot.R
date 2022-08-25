library("tidyr")
library("dplyr")
library("ggplot2")

technology <- "microarray"
data_set <- "GSE19188_m=21407"
path <- sprintf("results/diff-expr_%s_sample-size", technology)
fig_path <- "figures"

pattern <- sprintf("%s_SNR=(.*)_pi0=(.*)_N=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)

## check result availability
if (length(filenames) == 0L) {
  msg <- paste("Experiment results not found on disk. Please run", 
               "\t\tsource('scripts/sample-size-RNAseq_run.R')", 
               "\tto generate these results", 
               sep = "\n")
  stop(msg)
}

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
                  N = N, nb_exp = nb_exp)
  levList[[filename]] <- level
  
  power <- dat$power %>% group_by(alpha, method, selection) %>%
    summarise(estimate = mean(power, na.rm = TRUE),
              n = sum(!is.na(power)), # to check if large enough support for estimation
              std_err = sd(power, na.rm = TRUE)/sqrt(n),
              min = estimate - 2*std_err,
              max = estimate + 2*std_err,
              .groups = "drop")
  power <- tibble(power, pi0 = pi0, SNR = SNR, 
                  N = N, nb_exp = nb_exp)
  powList[[filename]] <- power
}

sc_pct <- function(x, ...) scales::percent(x, accuracy = 1, ...)

# - - - - - - - - - - - - - - - 
# Figure S9: JER contol
# - - - - - - - - - - - - - - - 
level <- Reduce(rbind, levList) 
dim(level)

level$method[level$method == "cherry/ARI"] <- "ARI"
level$method[level$method == "Simes (parametric)"] <- "Simes"
level$method[level$method == "Simes + step-down calibration"] <- "Adaptive Simes"
level$method[level$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

lev <- filter(level, 
              alpha <= 0.7, 
              pi0 %in% c(0.8), 
              N %in% c(10, 30, 50))
nb_exp_ <- unique(lev[["nb_exp"]])

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

plotname <- sprintf("fig-S9_JER-control_%s_%s_sample-size.pdf", 
                    technology, data_set)
pathname <- file.path(fig_path, plotname)
ggsave(p + geom_ribbon(alpha = 0.3, linetype = 1), 
       file = pathname, scale = 1, width = 16, height = 8)


# - - - - - - - - - - - - - - - 
# Power (not kept int he paper)
# - - - - - - - - - - - - - - - 
power <- Reduce(rbind, powList) 
dim(power)

power$method[power$method == "cherry/ARI"] <- "ARI"
power$method[power$method == "Simes (parametric)"] <- "Simes"
power$method[power$method == "Simes + step-down calibration"] <- "Adaptive Simes"
power$method[power$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

pow <- filter(power, 
              # estimate <= 1.0000001,
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

# plotname <- sprintf("fig-S9_power_%s_%s_sample-size.pdf", 
#                     technology, data_set)
# pathname <- file.path(fig_path, plotname)
# ggsave(p, file = pathname, scale = 1, width = 8, height = 6)
