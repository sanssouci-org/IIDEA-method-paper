library("sanssouci")
library("dplyr")
library("ggplot2")
library("DESeq2")
library("edgeR")
library("limma")

# - - - - - - - - - - 
# load RNAseq data
# - - - - - - - - - - 
data("RNAseq_blca", package = "sanssouci.data")
ds_name <- "BLCA"
X <- round(RNAseq_blca) # DESeq2 requires integers
groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
BLCA0 <- X/colSums(X)*1e6
ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
if (length(ww) != 0){
  X <- X[-ww, ]
}
rm(RNAseq_blca)
rm(BLCA0)

# - - - - - - - - - - 
# Wilcoxon tests
# - - - - - - - - - - 
wilcox <- rowWilcoxonTests(X, groups)

# - - - - - - - - - - 
# limma voom
# - - - - - - - - - - 
d <- DGEList(X)
# "Repeated column names found in count matrix"
d <- calcNormFactors(d)
Grp <- as.factor(groups)
mm <- model.matrix(~0 + Grp)
y <- voom(d, mm, plot = FALSE)
res_lm <- lmFit(y, mm)
contr <- makeContrasts(Grp1 - Grp0, levels = colnames(coef(res_lm)))
res_fit <- contrasts.fit(res_lm, contr)
res_eb <- eBayes(res_fit)
TT <- topTable(res_eb, sort.by = "none", number = Inf)
# plot(-log10(wilcox$p.value), -log10(TT$P.Value))

# - - - - - - - - - - 
# edgeR
# - - - - - - - - - - 
design <- model.matrix(~Grp)
d <- estimateDisp(d, design)
fit <- glmQLFit(d, design)
results <- glmQLFTest(fit)
edgeR <- topTags(results, sort.by = "none", n = Inf)
# plot(-log10(wilcox$p.value), -log10(edgeR$table$PValue))

# - - - - - - - - - - 
# DESeq
# - - - - - - - - - - 
dds <- DESeqDataSetFromMatrix(countData=X, colData = DataFrame(Grp),
                              design=~Grp, tidy = FALSE)
# "converting counts to integer mode"
dds
dds <- DESeq(dds)
deseq <- results(dds)
# plot(-log10(wilcox$p.value), -log10(deseq$pvalue))

# - - - - - - - - - - 
# gathering the results
# - - - - - - - - - - 
dat <- data.frame(wilcoxon = -log10(wilcox$p.value), 
                  limma_voom = -log10(TT$P.Value), 
                  edgeR = -log10(edgeR$table$PValue), 
                  DESeq = -log10(deseq$pvalue))
# - - - - - - - - - - 
# plots
# - - - - - - - - - - 
p1 <-  dat %>% 
  tidyr::pivot_longer(!wilcoxon, names_to = "test", values_to = "p_values") %>% 
  ggplot(aes(x = p_values, y = wilcoxon)) + 
  geom_point(alpha = 0.2, size = 0.5 ) + 
  facet_grid(.~test) + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "-log10(p-values)", y = "-log10(p_values) [wilcoxon]")
p1

p2 <- dat %>% 
  mutate(wilcoxon2 = wilcoxon,
         limma_voom2 = limma_voom, 
         edgeR2 = edgeR, 
         DESeq2 = DESeq) %>% 
  tidyr::pivot_longer(!c(wilcoxon, limma_voom, edgeR, DESeq), 
                      names_to = "test_x", 
                      values_to = "p_values_x") %>% 
  tidyr::pivot_longer(!c(test_x, p_values_x), 
                      names_to = "test_y", 
                      values_to = "p_values_y") %>%
  ggplot(aes(x = p_values_x, y = p_values_y)) + 
  geom_point(alpha = 0.2, size = 0.5 ) + 
  facet_grid(test_x~test_y) + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "-log10(p-values) [wilcoxon]", y = "-log10(p_values)")
p2
ggsave(filename = "comparison_p-value_all_tests.pdf", plot = p2)

## calibration
m <- dim(X)[1]
B <- 1000
alpha = 0.1
system.time({perm <- get_perm(X, groups, B, rowTestFUN = rowWilcoxonTests, alternative = "two.sided")})
perm_p_value <- perm$p.value
calibration <- calibrate(p0 = perm_p_value, m = nrow(X), alpha = alpha)
piv_stat <- calibration$piv_stat  
lambda <- stats::quantile(piv_stat, alpha, type = 1)
thr <- t_linear(lambda, 1:m, m)

sanssouci:::posthoc_bound(p.values = wilcox$p.value, 
                          thr = thr, 
                          lab = "Simes", 
                          all = TRUE, 
                          what = c("FDP","TP", "FP", "TDP" )
) %>%
  mutate(test = "wilcoxon") %>% 
  bind_rows(sanssouci:::posthoc_bound(p.values = TT$P.Value, 
                                      thr = thr, 
                                      lab = "Simes", 
                                      all = TRUE, 
                                      what = c("FDP","TP", "FP", "TDP" )
  ) %>%
    mutate(test = "limma-voom")
  )%>%
  bind_rows(sanssouci:::posthoc_bound(p.values = edgeR$table$PValue, 
                                      thr = thr, 
                                      lab = "Simes", 
                                      all = TRUE, 
                                      what = c("FDP","TP", "FP", "TDP" )
  ) %>%
    mutate(test = "edgeR")
  )%>%
  bind_rows(sanssouci:::posthoc_bound(p.values = deseq$pvalue, 
                                      thr = thr, 
                                      lab = "Simes", 
                                      all = TRUE, 
                                      what = c("FDP","TP", "FP", "TDP" )
  ) %>%
    mutate(test = "DESeq")
  ) -> curveTest 
saveRDS(curveTest, paste("curveTest",B,".rds"))
# curveTest <- readRDS("curveTest 1000 .rds")
p <- curveTest %>%  
  filter(stat %in% c("TP","FDP")) %>%
  ggplot(aes(x = x, y = bound, color= test))+
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~ stat, scales = "free_y") + 
  ggplot2::labs(x = "Number of top features selected", 
                y = paste("Post hoc confidence bounds (alpha =",alpha,")"), color = "Tests") +
  ggplot2::theme_bw() + 
  ggplot2::theme(strip.background = NULL, legend.position = "bottom") +
  xlim(c(0,1000))+
  ggplot2::scale_color_manual(values = scales::hue_pal()(4)) 
p
plotname <- sprintf("curveFDP_several_tests.pdf")
ggsave(p, file = plotname, scale = 1, width = 8, height = 6)


obj <- SansSouci(Y = X, groups = groups)
obj <- sanssouci::fit(obj, alpha = 0.1, B = 10, rowTestFun = rowWilcoxonTests)
plot(obj)






########### volcano plots ##############


obj <- SansSouci(Y = X, groups = groups)

# res_limma <- fit(obj, B = 0, family = "Simes", alpha = 0.05, rowTestFUN = rowLimmaVoom2)
# pvalLV <- pValues(res_limma)
# fcLV <- foldChanges(res_limma)
## limma voom
d <- DGEList(X)
#> Repeated column names found in count matrix
d <- calcNormFactors(d)
Grp <- as.factor(groups)
mm <- model.matrix(~0 + Grp)
y <- voom(d, mm, plot = FALSE)
res_lm <- lmFit(y, mm)
contr <- makeContrasts(Grp1 - Grp0, levels = colnames(coef(res_lm)))
res_fit <- contrasts.fit(res_lm, contr)
res_eb <- eBayes(res_fit)
TT <- topTable(res_eb, sort.by = "none", number = Inf)
pvalLV <- TT$P.Value
fcLV <- TT$logFC

volcanoP <- function(obj, pvalLV, fcLV, method = "Limma voom"){
  res_wilcoxon <- fit(obj, B = 1000, family = "Simes", alpha = 0.1, rowTestFUN = rowWilcoxonTests)
  
  pvalWcx <- pValues(res_wilcoxon)
  thrWcx <- thresholds(res_wilcoxon)
  
  p = 1e-3
  q = 1
  r = 0.5
  cex = c(0.2, 0.6)
  col = c("#33333333", "#FF0000", "#FF666633")
  pch = 19
  ylim = NULL
  bounds = TRUE
  # pval <- x; rm(x);
  if (p < 1 && q < 1) {
    warning("Filtering both on p-values and BH-adjusted p-values")
  }
  m <- length(pvalLV)
  
  ## sanity checks
  stopifnot(length(fcLV) == m)
  
  logpLV <- -log10(pvalLV)
  adjpLV <- p.adjust(pvalLV, method = "BH")  ## adjusted p-values
  y_sel <- which((adjpLV <= q) &           ## selected by q-value
                   (pvalLV <= p))        ##          or p-value
  y_thr <- Inf
  if (length(y_sel) > 0) {
    y_thr <- min(logpLV[y_sel])       ## threshold on the log(p-value) scale
  }
  
  ## gene selections
  sel1 <- which(logpLV >= y_thr & fcLV >= r)
  sel2 <- which(logpLV >= y_thr & fcLV <= -r)
  sel12 <- sort(union(sel1, sel2))
  
  ## post hoc bounds in selections
  n1 <- length(sel1)
  FP1 <- maxFP(pvalWcx[sel1], thr = thrWcx)
  TP1 <- n1 - FP1
  FDP1 <- round(FP1/max(n1, 1), 2)
  
  n2 <- length(sel2)
  FP2 <- maxFP(pvalWcx[sel2], thr = thrWcx)
  TP2 <- n2 - FP2
  FDP2 <- round(FP2/max(n2, 1), 2)
  
  n12 <- length(sel12)
  FP12 <- maxFP(pvalWcx[sel12], thr = thrWcx)
  TP12 <- n12 - FP12
  FDP12 <- round(FP12/max(n12, 1), 2)
  
  ## graphical parameters
  cols <- rep(col[1], m)
  cols[c(sel1, sel2)] <- col[2]
  
  cexs <- rep(cex[1], m)
  cexs[sel12] <- cex[2]
  
  xlab <- "Fold change (log scale)"
  ylab <- bquote("p-value (-" ~ log[10] ~ "scale)")
  infty <- 100
  if (is.null(ylim)) {
    ylim <- c(0, max(logpLV))
  }
  plot(fcLV, logpLV, pch = pch, cex = cexs, col = cols, 
       xlab = xlab, ylab = ylab, ylim = ylim)
  # axis4 <- thrYaxis(thr, max(ylim))
  # axis(side = 4, at = axis4$pvalue, labels = axis4$num, las = 1)
  # abline(h = -log10(thr[1:100]), col = "lightgray")
  rect(xleft = -infty, ybottom = y_thr, xright = -r, ytop = infty, 
       col = col[3], border = NA, lwd = 2)
  rect(xleft = r, ybottom = y_thr, xright = infty, ytop = infty, 
       col = col[3], border = NA, lwd = 2)
  abline(h = y_thr, col = "gray")
  abline(v = c(-1, 1)*r, col = "gray")
  
  if (bounds) {
    bq <- bquote(atop(.(n1) ~ "genes", 
                      "TP"  >= .(TP1) ~ ";  FDP" <= .(FDP1)))
    legend("topright", legend = bq, border = "white", bty = "n", text.col = 1)
    bq <- bquote(atop(.(n2) ~ "genes",
                      "TP"  >= .(TP2) ~ "; FDP" <= .(FDP2)))
    legend("topleft", legend = bq, border = "white", bty = "n", text.col = 1)
    
    bq <- bquote(atop(.(method) ~ ":" ~ .(n12) ~ "genes selected",
                      "At least" ~ .(TP12) ~ "true positives (FDP" <= .(FDP12) ~")"))
    title(bq)
  }
  invisible(sel12)
}
volcanoP(obj, pvalLV = pvalLV, fcLV = fcLV, method = "Limma voom")

design <- model.matrix(~Grp)
d <- estimateDisp(d, design)
fit <- glmQLFit(d, design)
results <- glmQLFTest(fit)
edgeR <- topTags(results, sort.by = "none", n = Inf)
volcanoP(obj, pvalLV = edgeR$table$PValue, fcLV = edgeR$table$logFC, method = "EdgeR")

## DESeq
dds <- DESeqDataSetFromMatrix(countData=X, colData = DataFrame(Grp),
                              design=~Grp, tidy = FALSE)
dds
dds <- DESeq(dds)
deseq <- results(dds)
volcanoP(obj, pvalLV = deseq$pvalue, fcLV = deseq$log2FoldChange, method = "DESeq2")

par(mfrow = c(1,3))
volcanoP(obj, pvalLV = pvalLV, fcLV = fcLV, method = "Limma voom")
volcanoP(obj, pvalLV = edgeR$table$PValue, fcLV = edgeR$table$logFC, method = "EdgeR")
volcanoP(obj, pvalLV = deseq$pvalue, fcLV = deseq$log2FoldChange, method = "DEseq")

ggplot(data = data.frame(fc = fcLV, logpvalues = logpLV, cex = cexs, col = cols)) + 
  geom_point(aes(x = fc, y = logpvalues, size = cex, color = col)) + 
  geom_vline(xintercept = c(-1, 1)*r, col = "gray") + 
  geom_hline(yintercept = y_thr, col = "gray") +
  annotate("rect", xmin = -min(fcLV), xmax = -r, ymin = y_thr, ymax = infty,
           alpha = .1,fill = col[3])  + 
  ylim(c(0, max(logpLV)))  + 
  xlim(c(min(fcLV), max(fcLV))) 

res_i$power %>% select(time) %>% unique()


