geo2kegg <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData, "geo2kegg")
idx <- match(ds_name, names(geo2kegg))
if (is.na(idx)) {
  stop("Data set not found in geo2kegg: ", ds_name)
}
ds <- geo2kegg[idx]
rawData <- R.cache::memoizedCall(maPreproc, ds)[[1]]
X <- SummarizedExperiment::assays(rawData)$exprs
cats <- SummarizedExperiment::colData(rawData)
ww <- match(cats$Sample, base::colnames(X))
groups <- cats$GROUP[ww]

X0 <- X[, groups == 0]
