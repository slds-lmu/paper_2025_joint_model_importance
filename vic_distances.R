# Pairwise importance distances between Rashomon-set models.
#
# Mirrors main.R's behavioural-distance pipeline but uses the per-model
# feature-importance vectors stored in data/results_vic.RData (the VIC,
# Variable Importance Cloud) instead of the validation predictions. The
# point is to be able to compare the two views: behavioural distance
# answers "do these models predict the same?", VIC distance answers
# "do these models rank features the same?". The two are independent
# (identical predictions need not imply identical importance) so a
# behavioural cluster need not coincide with a VIC cluster -- that
# coincidence (or its absence) is the joint-importance question.
#
# Reads:
#   data/run_models_merged.rds (via get_performance_and_SVMkernel)
#   data/results_vic.RData     (loaded by init/zzz_settings.R as `vic`)
#   data/preds_cache.rds       (optional, for the side-by-side plot)
#
# Writes (figures/vic/):
#   dendrogram_{task}_{metric}.pdf
#   mds_{task}_{metric}.pdf
#   pam_{task}_{metric}.pdf
#   mds_compare_{task}.pdf      (behavioural vs VIC, Euclidean only)

source("init.R")

run_models_merged_location <- "data/run_models_merged.rds"
preds_cache_location       <- "data/preds_cache.rds"
figures_dir                <- "figures"
vic_dir                    <- file.path(figures_dir, "vic")

RS_epsilon       <- 0.01
distance_metrics <- c("euclidean", "manhattan")

dir.create(vic_dir, showWarnings = FALSE, recursive = TRUE)

# Rashomon-set membership per task (same call main.R uses)
perfor_and_kernel <- get_performance_and_SVMkernel(run_models_merged_location)
performance       <- perfor_and_kernel$performance
RS <- get_RS(RS_epsilon, performance, vic)

# Optional: behavioural prediction cache, used for the comparison plot
preds_by_task <- if (file.exists(preds_cache_location)) {
  readRDS(preds_cache_location)
} else {
  message("No preds cache at ", preds_cache_location,
          " -- side-by-side comparison plots will be skipped.")
  NULL
}

# VIC column names follow "pfi_{learner}_m{model_no}" (see init/functions.R)
parse_vic_colnames <- function(cnames) {
  list(
    learner  = sub(".*_(.*?)_.*", "\\1", cnames),
    model_no = as.integer(sub(".*_m([0-9]+)$", "\\1", cnames))
  )
}

# Per task: importance matrix with one row per RS model, one column per
# feature; plus the aligned learner / model_no vectors.
# vic[[task]] has the feature names in column 1; RS[[task]] gives row
# indices in design order, which match the column order in vic past the
# feature column. Hence the +1 offset.
vic_by_task <- lapply(names(RS), function(task_name) {
  full <- vic[[task_name]]
  if (is.null(full)) return(NULL)
  rs_cols <- RS[[task_name]] + 1L
  imp_mat <- t(as.matrix(full[, rs_cols, drop = FALSE]))
  parsed  <- parse_vic_colnames(colnames(full)[rs_cols])
  rownames(imp_mat) <- sprintf("%s_m%d", parsed$learner, parsed$model_no)
  list(
    imp      = imp_mat,
    learner  = parsed$learner,
    model_no = parsed$model_no,
    features = full[[1]]
  )
})
names(vic_by_task) <- names(RS)

# Sanity check: VIC and behavioural model orderings should agree so the
# comparison plot colours line up correctly.
if (!is.null(preds_by_task)) {
  for (task_name in names(vic_by_task)) {
    td <- vic_by_task[[task_name]]; pb <- preds_by_task[[task_name]]
    if (is.null(td) || is.null(pb)) next
    if (length(td$model_no) != length(pb$model_no) ||
        any(td$model_no != pb$model_no)) {
      warning("VIC and preds model orderings differ for task ", task_name,
              " -- comparison plot may be misaligned.")
    }
  }
}

#### per-task plots ###########################################################

for (task_name in names(vic_by_task)) {
  td <- vic_by_task[[task_name]]
  if (is.null(td) || nrow(td$imp) < 2) next
  learner_factor <- factor(td$learner)

  for (metric in distance_metrics) {
    dist_mat <- dist(td$imp, method = metric)

    # dendrogram
    pdf(file.path(vic_dir, sprintf("dendrogram_%s_%s.pdf", task_name, metric)),
        width = 9, height = 6)
    plot(hclust(dist_mat), hang = -1,
         main = sprintf("VIC dendrogram -- %s (%s)", task_name, metric),
         xlab = "model", sub = "", cex = 0.6)
    dev.off()

    # MDS coloured by learner
    mds <- cmdscale(dist_mat, k = 2, eig = TRUE)
    pdf(file.path(vic_dir, sprintf("mds_%s_%s.pdf", task_name, metric)),
        width = 7, height = 6)
    plot(mds$points, type = "n",
         xlab = "MDS 1", ylab = "MDS 2", asp = 1,
         main = sprintf("VIC MDS -- %s (%s)", task_name, metric))
    points(mds$points,
           col = as.integer(learner_factor), pch = 19, cex = 0.9)
    legend("topright",
           legend = levels(learner_factor),
           col    = seq_len(nlevels(learner_factor)),
           pch    = 19, bty = "n")
    dev.off()

    # PAM placeholder (same k=2 convention as the behavioural pipeline)
    pam_res <- cluster::pam(dist_mat, k = 2, metric = metric)
    pdf(file.path(vic_dir, sprintf("pam_%s_%s.pdf", task_name, metric)),
        width = 9, height = 5)
    plot(pam_res, main = sprintf("VIC PAM -- %s (%s)", task_name, metric))
    dev.off()
  }
}

#### side-by-side comparison (Euclidean only) #################################
# Behavioural MDS on the left, VIC MDS on the right, same task, same colour
# scheme (one colour per learner family). Whether the two views agree is the
# joint-importance question.

if (!is.null(preds_by_task)) {
  for (task_name in names(vic_by_task)) {
    td <- vic_by_task[[task_name]]; pb <- preds_by_task[[task_name]]
    if (is.null(td) || is.null(pb)) next

    bh_mds <- cmdscale(dist(pb$preds, method = "euclidean"), k = 2)
    vc_mds <- cmdscale(dist(td$imp,   method = "euclidean"), k = 2)
    bh_lf  <- factor(pb$learner)
    vc_lf  <- factor(td$learner)

    pdf(file.path(vic_dir, sprintf("mds_compare_%s.pdf", task_name)),
        width = 13, height = 6)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(bh_mds, col = as.integer(bh_lf), pch = 19, cex = 0.9, asp = 1,
         xlab = "MDS 1", ylab = "MDS 2",
         main = sprintf("Behavioural distance -- %s", task_name))
    legend("topright", legend = levels(bh_lf),
           col = seq_len(nlevels(bh_lf)), pch = 19, bty = "n", cex = 0.8)
    plot(vc_mds, col = as.integer(vc_lf), pch = 19, cex = 0.9, asp = 1,
         xlab = "MDS 1", ylab = "MDS 2",
         main = sprintf("VIC distance -- %s", task_name))
    legend("topright", legend = levels(vc_lf),
           col = seq_len(nlevels(vc_lf)), pch = 19, bty = "n", cex = 0.8)
    dev.off()
  }
}

message("Wrote VIC plots to ", vic_dir, "/")
