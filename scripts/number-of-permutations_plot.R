library("tidyr")
library("dplyr")
library("ggplot2")

technology <- "RNAseq"
# data_set <- "BLCA_m=12418"
data_set <- "BLCA_m=1234"
path <- sprintf("results/diff-expr_%s_permutation", technology)
fig_path <- "figures"

pattern <- sprintf("%s_B=(.*)_nb-exp=(.*).rds",
                   data_set)
filenames <- list.files(path, pattern = pattern)

## check result availability
if (length(filenames) == 0L) {
    msg <- paste("Experiment results not found on disk. Please run", 
                 "\t\tsource('scripts/number-of-permutations.R')", 
                 "\tto generate these results", 
                 sep = "\n")
    stop(msg)
}

lambdaList <- list()
FDPList <- list()
for (filename in filenames) {
    print(filename)
    pathname <- file.path(path, filename)
    B <- as.numeric(gsub(pattern, "\\1", filename))
    nb_exp <- as.numeric(gsub(pattern, "\\2", filename))
    dat <- readRDS(pathname)
    
    lambda <- dat$lambda %>% group_by(alpha, method) 
    lambda <- tibble(lambda, B = B, nb_exp = nb_exp)
    lambdaList[[filename]] <- lambda
    
    FDP <- dat$FDP %>% group_by(alpha, method, selection) 
    FDP <- tibble(FDP, B = B, nb_exp = nb_exp)
    FDPList[[filename]] <- FDP
}

## FDP
FDP <- Reduce(rbind, FDPList) %>%
    filter(alpha == 0.1, 
           selection %in% c("BH_05", "first_1000", "H")
    ) %>%
    mutate(B = factor(B)) %>%
    rename(FDP = power)


p <- FDP %>% 
    filter(selection %in% c("BH_05", "first_1000", "H")) %>%
    group_by(alpha, selection, B) %>% 
    summarise(mean_time = mean(time), 
              q_01 = quantile(FDP, probs = 0.01), 
              q_99 = quantile(FDP, probs = 0.99)) %>% 
    mutate(diff_q = q_99 - q_01) %>%
    filter(alpha == 0.1) %>%
    ggplot(aes(x = mean_time, y = diff_q, color = selection, label = B)) + 
    geom_point() + 
    geom_line() +
    geom_text(hjust=-0.1, vjust=-0.1) +
    theme_minimal() +
    xlab("Mean time (across 1000 experiments)") + ylab("Inter-centile range (Q0.99 - Q0.01) (across 1000 experiments)")
p

plotname <- sprintf("fig-S4_permutation-time-vs-ICR_%s_%s.pdf", technology, data_set)
print(plotname)
pathname <- file.path(fig_path, plotname)
ggsave(p, file = pathname, scale = 1, width = 8, height = 6)

B <- 1000

FDP %>% 
    filter(selection %in% c("BH_05", "first_1000", "H")) %>%
    group_by(alpha, selection, B) %>% 
    summarise(mean_time = mean(time)) %>% filter(B == B) %>% ungroup() %>% select(mean_time) %>% unique()

FDP %>% filter(B == 100, alpha == 0.1, selection == "BH_05") %>%
    group_by(B) %>% 
    summarise(q0.995 = quantile(FDP, probs = 0.995), 
              q0.005 = quantile(FDP, probs = 0.005), 
              q0.99 = quantile(FDP, probs = 0.99),
              q0.01 = quantile(FDP, probs = 0.01))

p <- FDP %>%
    filter(alpha == 0.1, 
           selection %in% c("BH_05", "first_1000", "H")
    ) %>%
    # mutate(B = factor(B)) %>% 
    mutate(FDP = FDP)%>% 
    left_join(FDP %>% group_by(B) %>% summarise(t = mean(time)), by="B")%>%
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

plotname <- "fig-S3_FDP-distribution_permutations_BLCA_m=12534.pdf"
print(plotname) 
pathname <- file.path(fig_path, plotname)
ggsave(p, file = pathname, scale = 1, width = 8, height = 6)
