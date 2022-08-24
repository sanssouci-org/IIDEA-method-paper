data("RNAseq_blca", package = "sanssouci.data")
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

rm(X)
