#!/usr/bin/env Rscript
# R reference imputers for the v9 addendum (post hoc, exploratory). Called by src/scp_graphfs_mcp/addendum2/rbridge.py.
#
# Usage: Rscript scripts/r_bridge.R <in.csv|in.rds> <method> <seed> <out.csv>
#   in   cells x proteins matrix (log2), NA = missing; CSV without header/row names, or an RDS numeric matrix.
#   out  the imputed cells x proteins matrix, same shape, CSV with 17 significant digits, no header.
#   The caller must drop proteins (columns) and cells (rows) that are entirely NA (pcaMethods::checkData refuses them).
#
# Methods (orientation follows the reference wrappers in MsCoreUtils::impute_matrix, used by MSnbase::impute, where the
# matrix is features x samples, i.e. proteins x cells):
#   bpca        pcaMethods::pca(t(X), method = "bpca", nPcs = ncol(t(X)) - 1, verbose = FALSE)  (MsCoreUtils::impute_bpca)
#   lls         pcaMethods::llsImpute(X, k = 10, completeObs = TRUE)  - pcaMethods wants variables (genes) in columns, so
#               cells x proteins; every other argument at its default (it switches to allVariables = TRUE itself when
#               fewer than k proteins are complete, which is always the case here)
#   lls_onepass NOT a default: llsImpute(X, k = 10, center = TRUE, maxSteps = 1). With defaults (center = FALSE,
#               maxSteps = 100, unregularised ginv on ~20 observed cells x k = 10 neighbours) the allVariables iteration
#               diverges on 30-cell matrices (synthetic RMSE ~1e10); this variant stops after the first pass.
#   knn         impute::impute.knn(t(X), k = 10, rng.seed = seed)$data, other arguments default (MsCoreUtils::impute_knn)
#   qrilc       imputeLCMD::impute.QRILC(t(X), tune.sigma = 1)[[1]]                          (MsCoreUtils::impute_QRILC)
#   missforest  missForest::missForest(t(X), num.threads = 1)$ximp: proteins are the observations and cells the variables
#               (MsCoreUtils::impute_RF orientation); defaults otherwise (ntree 100, maxiter 10, ranger backend)
#   missforest_cells  missForest::missForest(X, num.threads = 1)$ximp: cells are the observations, proteins the variables
#   msimpute    msImpute::msImpute(t(X), method = "v2-mnar", group = rep(1, ncol(t(X))))  - the benchmark has no group
#               labels at imputation time, so one group (an assumption, not a msImpute default; group is required)
#   msimpute_v2 msImpute::msImpute(t(X), method = "v2")   (MAR-only low-rank part, no group needed)

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4) stop("usage: r_bridge.R <in.csv|in.rds> <method> <seed> <out.csv>")
})
inp <- args[1]; method <- args[2]; seed <- as.integer(args[3]); outp <- args[4]

if (grepl("\\.rds$", inp, ignore.case = TRUE)) {
  X <- readRDS(inp)
} else {
  X <- as.matrix(utils::read.csv(inp, header = FALSE, na.strings = c("NA", "nan", "NaN", "")))
}
storage.mode(X) <- "double"
X[is.nan(X)] <- NA
dimnames(X) <- list(paste0("c", seq_len(nrow(X))), paste0("p", seq_len(ncol(X))))
set.seed(seed)

P <- t(X)  # proteins x cells (features x samples)
Y <- switch(method,
  bpca = {
    res <- pcaMethods::pca(P, method = "bpca", nPcs = ncol(P) - 1, verbose = FALSE)
    t(pcaMethods::completeObs(res))
  },
  lls = {
    res <- suppressWarnings(pcaMethods::llsImpute(X, k = 10, completeObs = TRUE))
    pcaMethods::completeObs(res)
  },
  lls_onepass = {
    res <- suppressWarnings(pcaMethods::llsImpute(X, k = 10, center = TRUE, maxSteps = 1, completeObs = TRUE))
    pcaMethods::completeObs(res)
  },
  knn = {
    t(impute::impute.knn(P, k = 10, rng.seed = seed)$data)
  },
  qrilc = {
    t(imputeLCMD::impute.QRILC(P, tune.sigma = 1)[[1]])
  },
  missforest = {
    t(missForest::missForest(P, num.threads = 1)$ximp)
  },
  missforest_cells = {
    missForest::missForest(X, num.threads = 1)$ximp
  },
  msimpute = {
    suppressMessages(t(msImpute::msImpute(P, method = "v2-mnar", group = rep(1, ncol(P)))))
  },
  msimpute_v2 = {
    suppressMessages(t(msImpute::msImpute(P, method = "v2")))
  },
  stop(paste("unknown method", method))
)
Y <- as.matrix(Y)
if (!all(dim(Y) == dim(X))) stop("output shape mismatch")
Y <- Y[rownames(X), colnames(X), drop = FALSE]
lines <- apply(Y, 1, function(r) paste(ifelse(is.finite(r), sprintf("%.17g", r), "NA"), collapse = ","))
writeLines(lines, outp)
