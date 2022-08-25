library("tidyr")
library("dplyr")
library("ggplot2")

technology <- "simulation"
data_set <- "block-cov_m=1000"
path <- file.path("results", sprintf("diff-expr_%s_correlation", technology))
pattern <- sprintf("%s_nb_cor=(.*)_rho=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)

## check result availability
if (length(filenames) == 0L) {
    msg <- paste("Experiment results not found on disk. Please run", 
                 "\t\tsource('scripts/continuous-covariate_run.R')", 
                 "\tto generate these results", 
                 sep = "\n")
    stop(msg)
}

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
lev <- filter(level, alpha <=0.7, 
              nb_cor %in% c(1, 5, 20),
              rho %in% c(0.2, 0.4, 0.6))
nb_exp_ <- unique(lev[["nb_exp"]])

sc_pct <- function(x, ...) scales::percent(x, accuracy = 1, ...)

# - - - - - - - - - - - - - - - -
# Figure S-10: JER control
# - - - - - - - - - - - - - - - -
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
plotname <- sprintf("fig-S10_JER-control_%s_%s_correlation.pdf", technology, data_set)
ggsave(p + geom_ribbon(alpha = 0.3, linetype = 1), 
       file = plotname, scale = 1, width = 8, height = 8)

# - - - - - - - - - - - - - - - -
# Power (not shown in the paper)
# - - - - - - - - - - - - - - - -

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
# p
# plotname <- sprintf("power_%s_%s_correlation.pdf", technology, data_set)
# ggsave(p, file = plotname, scale = 1, width = 8, height = 6)
