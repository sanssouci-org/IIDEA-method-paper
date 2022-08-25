# - - - - - - - - - - - - -
# setup
# - - - - - - - - - - - - -
library("ggplot2")
library("tidyr")
library("sanssouci")
library("matrixStats")
fig_path <- "figures"
dir.create(fig_path, showWarnings = FALSE)

# - - - - - - - - - - - - -
# data
# - - - - - - - - - - - - -
data("RNAseq_blca", package = "sanssouci.data")
Y <- RNAseq_blca
groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
rm(RNAseq_blca)
dim(Y)

CPM <- Y/colSums(Y)*1e6
ww <- which(rowQuantiles(CPM, prob = 0.75) < 5)
if (length(ww) != 0){
    Y <- Y[-ww, ]
}
rm(CPM, ww)
dim(Y)

set.seed(19012001)

alpha <- 0.1

# - - - - - - - - - - - - - - -
# Adaptive Simes (step down)
# - - - - - - - - - - - - - - -
obj <- SansSouci(Y = log(1 + Y), groups = groups)
res <- fit(obj, B = 1000, alpha = alpha, family = "Simes", 
           rowTestFUN = rowWilcoxonTests)

# - - - - - - - - - - - - - - -
# Adaptive Simes (single step)
# - - - - - - - - - - - - - - -
res_singlestep <- fit(obj, B = 1000, alpha = alpha, family = "Simes", 
                      rowTestFUN = rowWilcoxonTests, max_steps_down = 0)

# - - - - - - - - - - - - - - -
# Parametric Simes (single step)
# - - - - - - - - - - - - - - -
res_Simes <- fit(obj, B = 0, family = "Simes", alpha = alpha, 
                 rowTestFUN = rowWilcoxonTests) ## B=0 => no calibration!

m <- dim(Y)[1]
d <- hommel::discoveries(hommel::hommel(pValues(res)))
pi0_hat <- 1 - d/m

# - - - - - - - - - - - - - - -
# Parametric Simes (step down) aka ARI
# - - - - - - - - - - - - - - -
res_ARI <- fit(obj, B = 0, family = "Simes", alpha = alpha/pi0_hat, 
               rowTestFUN = rowWilcoxonTests)

resList <- list("Adaptive Simes" = res,
                "Adaptive Simes (single step)" = res_singlestep,
                "Simes" = res_Simes,
                "ARI" = res_ARI)

# - - - - - - - - - - - - - - -
# Calculating posthoc bounds
# - - - - - - - - - - - - - - -
m <- nrow(Y)
q <- 0.1 # FDP budget (user-defined)
FDP <- lapply(resList, predict, what = "FDP", all = TRUE)
n_DEG <- sapply(FDP, function(x) max(which(x$bound <= q)))

TP <- NULL
for (method in names(resList)) {
    n <- n_DEG[[method]]
    TP[method] <- predict(resList[[method]], what = "TP", all = TRUE)$bound[n]
}

DEG_all <- tibble(Method = names(FDP), x = n_DEG, TP = TP, FDP = q) %>%
    pivot_longer(c(TP, FDP), names_to = "stat", values_to = "y")
conf_bounds_all <- lapply(resList, predict, all = TRUE)


# - - - - - - - - - - - - - - -
# Genes selected by BH
# - - - - - - - - - - - - - - -
bh <- which(p.adjust(pValues(res), method = "BH") < 0.05)
xBH <- length(bh)
tpBH_res <- predict(res, what = c("TP"), all = TRUE)$bound[xBH]
fdpBH_res <- predict(res, what = c("FDP"), all = TRUE)$bound[xBH]
tpBH_resSimes <- predict(res_ARI, what = c("TP"), all = TRUE)$bound[xBH]
fdpBH_resSimes <- predict(res_ARI, what = c("FDP"), all = TRUE)$bound[xBH]

BH <- tibble(Template = "Adaptive Simes", 
             x = xBH, 
             TP = c(tpBH_res), 
             FDP = c(fdpBH_res)) %>%
    pivot_longer(c(TP, FDP), names_to = "stat", values_to = "y")

# - - - - - - - - - - - - - - -
# Figure 2 (confidence curves)
# - - - - - - - - - - - - - - -
methods <- c("Adaptive Simes", "ARI")
DEG <- subset(DEG_all, Method %in% methods)
conf_bounds <- conf_bounds_all[methods]

cols <- c("black", "lightgray")
p <- plotConfCurve(conf_bounds, xmax = 2.5*max(n_DEG), cols = cols) +  
    scale_linetype_manual(values = c("solid", "solid")) +
    labs(color = "Method", linetype = "Method") +
    #    geom_vline(xintercept = n_DEG, linetype = "dotted", size = size) + 
    geom_segment(data = DEG, aes(x = x, y = -Inf, xend = x, yend = y, 
                                 color = Method, linetype = Method), 
                 size = 1) + 
    geom_segment(data = DEG, aes(x = -Inf, y = y, xend = x, yend = y, 
                                 color = Method, linetype = Method), 
                 size = 1) +
    geom_line(size = 1.2) + 
    geom_point(data = BH, aes(x = x, y = y), colour = "red", size = 2.5) +
    ggplot2::facet_wrap(~ stat, scales = "free_y") + 
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = xBH, linetype="dotted", 
               color = "red", size=0.5) 
#   + geom_hline(yintercept = q, linetype = "dashed", size = size) 
p

plotname <- "fig-2_conf-curve.pdf"
pathname <- file.path(fig_path, plotname)
ggsave(p, file = pathname, width = 6, height = 4)

# - - - - - - - - - - - - - - -
# Figure S-1 (all methods)
# - - - - - - - - - - - - - - -
conf_bounds <- conf_bounds_all
DEG <- DEG_all

cols <- rep(c("black", "lightgray"), each = 2)
ltys <- rep(c("solid", "dashed"), times = 2)

p <- plotConfCurve(conf_bounds, xmax = 2.5*max(n_DEG), cols = cols) +  
    scale_linetype_manual(values = ltys) +
    labs(color = "Method", linetype = "Method") +
    geom_line(size = 1.2) + 
    ggplot2::facet_wrap(~ stat, scales = "free_y") +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) 
p

plotname <- "fig-S1_conf-curve-annexe.pdf"
pathname <- file.path(fig_path, plotname)
ggsave(p, file = pathname, width = 6, height = 4)

# - - - - - - - - - - - - - - -
# Figure 3 (volcano plot)
# - - - - - - - - - - - - - - -
m <- nrow(Y)
q <- 0.1 # FDP budget (user-defined)
FDP <- lapply(resList, predict, what = "FDP", all = TRUE)
n_DEG <- sapply(FDP, function(x) max(which(x$bound <= q)))

library(limma)
library(edgeR)
d <- DGEList(Y)
d <- calcNormFactors(d)
Grp <- as.factor(groups)
mm <- model.matrix(~0 + Grp)

y <- voom(d, mm, plot = FALSE)

res_lm <- lmFit(y, mm)
contr <- makeContrasts(Grp1 - Grp0, levels = colnames(coef(res_lm)))
res_fit <- contrasts.fit(res_lm, contr)
res_eb <- eBayes(res_fit)
TT <- topTable(res_eb, sort.by = "none", number = Inf)

plotname <- "fig-3_volcano-plot.pdf"
pathname <- file.path(fig_path, plotname)
pdf(pathname, width = 6, height = 6)
volcanoPlot(res, 
            fold_changes = TT$logFC, 
            p_values = TT$P.Value, 
            p = 1e-3, r = 0.5)
dev.off()


# - - - - - - - - - - - - - - - - - - - -
# Table 2 (right part)
# - - - - - - - - - - - - - - - - - - - -
selVP <- volcanoPlot(res, 
                     fold_changes = TT$logFC, 
                     p_values = TT$P.Value, 
                     p = 1e-3, r = 0.5)
TP <- list()
for (method in names(resList)) {
    n <- n_DEG[[method]]
    TP[method] <- list(predict(resList[[method]], what = c("TP", "FDP"), S = selVP))
    #    c(Template = method, n = n, TP = TP)
}
length(selVP)
TP

# - - - - - - - - - - - - - - - - - - - -
# Figure S-2 (limma vs Wilcoxon p-values)
# - - - - - - - - - - - - - - - - - - - -
df <- data.frame(wilcox = -log10(pValues(res)), limma = -log10(TT$P.Value))
p <- ggplot(df, aes(x = limma, y = wilcox)) + 
    geom_point(color = "#10101010") + 
    ggtitle("p-values (log-scale)") + 
    theme_bw() + theme(legend.position="bottom")
p
plotname <- "fig-S2_p-values_limma-vs-wilcoxon.pdf"
pathname <- file.path(fig_path, plotname)
ggsave(p, file = pathname, width = 6, height = 4)
