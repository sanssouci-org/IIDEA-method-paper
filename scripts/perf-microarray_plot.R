library("tidyr")
library("dplyr")
library("ggplot2")

technology <- "microarray"
data_set <- "GSE19188_m=21407"
path <- sprintf("results/diff-expr_%s_perf", technology)
fig_path <- "figures"

pattern <- sprintf("%s_pi0=(.*)_SNR=(.*)_prob=(.*)_B=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)

## check result availability
if (length(filenames) == 0L) {
    msg <- paste("Experiment results not found on disk. Please run", 
                 "\t\tsource('scripts/perf-microarray_run.R')", 
                 "\tto generate these results", 
                 sep = "\n")
    stop(msg)
}

levList <- list()
powList <- list()
for (filename in filenames) {
    print(filename)
    pathname <- file.path(path, filename)
    pi0 <- as.numeric(gsub(pattern, "\\1", filename))
    SNR <- gsub(pattern, "\\2", filename)
    prob <- as.numeric(gsub(pattern, "\\3", filename))
    B <- as.numeric(gsub(pattern, "\\4", filename))
    nb_exp <- as.numeric(gsub(pattern, "\\5", filename))
    dat <- readRDS(pathname)
    
    level <- dat$level %>% group_by(alpha, method) %>%
        summarise(estimate = 1 - mean(`valid bound`), 
                  std_err = sqrt((1 - estimate)*estimate/n()),
                  min = estimate - 2*std_err,
                  max = estimate + 2*std_err,
                  .groups = "drop")
    level <- tibble(level, pi0 = pi0, SNR = SNR, prob = prob,
                    B = B, nb_exp = nb_exp)
    levList[[filename]] <- level
    
    power <- dat$power %>% group_by(alpha, method, selection) %>%
        summarise(estimate = mean(power, na.rm = TRUE),
                  n = sum(!is.na(power)), # to check if large enough support for estimation
                  std_err = sd(power, na.rm = TRUE)/sqrt(n),
                  min = estimate - 2*std_err,
                  max = estimate + 2*std_err,
                  .groups = "drop")
    power <- tibble(power, pi0 = pi0, SNR = SNR, 
                    prob = prob, B = B, nb_exp = nb_exp)
    powList[[filename]] <- power
}
level <- Reduce(rbind, levList) 
dim(level)

level$method[level$method == "cherry/ARI"] <- "ARI"
level$method[level$method == "Simes (parametric)"] <- "Simes"
level$method[level$method == "Simes + step-down calibration"] <- "Adaptive Simes"
level$method[level$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

lev <- filter(level, alpha <= 0.5, pi0 %in% c(0.5, 0.8, 0.9, 1), prob == 0.5)
nb_exp_ <- unique(lev[["nb_exp"]])

sc_pct <- function(x, ...) scales::percent(x, accuracy = 1, ...)

# - - - - - - - - - - - - - - - 
# Figure S-8: JER contol
# - - - - - - - - - - - - - - - 
p <- ggplot(lev, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method,
                     fill = method,
                     shape = method)) + 
    facet_grid(SNR~pi0, labeller = label_bquote(cols = pi[0]: .(pi0),
                                                rows = SNR: .(SNR),)) + 
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

plotname <- sprintf("fig-S6_JER-control_%s_%s.pdf", technology, data_set)
print(plotname)
pathname <- file.path(fig_path, plotname)
ggsave(p + geom_ribbon(alpha = 0.3, linetype = 1), 
       file = pathname, scale = 1, width = 8, height = 8)

# - - - - - - - - - - - - - - - 
# Power (not shown in the paper)
# - - - - - - - - - - - - - - - 
power <- Reduce(rbind, powList) 
dim(power)

power$method[power$method == "cherry/ARI"] <- "ARI"
power$method[power$method == "Simes (parametric)"] <- "Simes"
power$method[power$method == "Simes + step-down calibration"] <- "Adaptive Simes"
power$method[power$method == "Simes + single-step calibration"] <- "Adaptive Simes (single step)"

pow <- filter(power, 
              prob == 0.5,
              SNR %in% c("+2"),
              selection %in% c("first_100", "BH_05", "p_05", "H"), 
              alpha <= 0.5, 
              pi0 == 0.8)
nb_exp_ <- unique(pow[["nb_exp"]])
pow

p <- ggplot(pow, aes(x = alpha, y = estimate, 
                     ymax = max, ymin = min,
                     color = method, shape = method)) + 
    facet_wrap(~selection)+
    geom_point() + 
    geom_line() + 
    ylim(c(0, 1.0001)) +
    ylab(paste("Power (estimated from", nb_exp_, "experiments)")) +
    xlab(expression(paste("Target risk (", alpha, ")"))) +
    labs(color = "Method", fill = "Method", shape = "Method")+
    scale_x_continuous(labels = sc_pct) + 
    scale_y_continuous(labels = sc_pct) + 
    theme(legend.position="bottom")
# p
# plotname <- sprintf("power_%s_%s.pdf", technology, data_set)
# print(plotname)
# pathname <- file.path(fig_path, plotname)
# ggsave(p, file = pathname, scale = 1, width = 8, height = 6)
