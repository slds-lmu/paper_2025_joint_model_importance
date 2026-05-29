# Hyperparameter trajectories along the within-family MDS axis.
#
# Reads data/preds_cache.rds (written by main.R) and, for each
# (task, learner-family) with enough Rashomon-set members, loads the
# corresponding saved model files to read out the tuning hyperparameters.
# Then projects each model onto a 1D within-family MDS coordinate and plots
# that coordinate against every hyperparameter that varies across the family.
# One multi-panel PDF per (task, learner) into figures/hp_trajectories/.
#
# Model loading is slow (one readRDS per Rashomon-set member, hundreds of
# RDS reads in total), so the extracted HPs are cached to data/hp_cache.rds.
# Delete that file to force a reload.

source("init.R")
library(data.table)

preds_cache_location <- "data/preds_cache.rds"
hp_cache_location    <- "data/hp_cache.rds"
figures_dir          <- "figures"
hp_dir               <- file.path(figures_dir, "hp_trajectories")
model_root           <- "/media/external/rashomon/rashomon_models"
min_family_size      <- 5

if (!file.exists(preds_cache_location)) {
  stop("No predictions cache at ", preds_cache_location,
       " -- run main.R first.")
}
preds_by_task <- readRDS(preds_cache_location)

dir.create(hp_dir, showWarnings = FALSE, recursive = TRUE)

# Keys that show up in param_set$values but are bookkeeping rather than
# tuned hyperparameters; dropped from plots even if they happen to vary.
reserved_keys <- c("predict_type", "fitted", "type", "MaxNWts", "maxit",
                   "trace", "tolerance")

# For most learners (xgb, glmnet, nnet, svm) the saved model is a
# GraphLearner whose tuned HPs are namespaced "<learner>.<hp>". Pull only
# those and drop the prefix. tree models are not graph-wrapped, so HPs
# arrive without prefix -- fall through to using all of them.
extract_hps_from_model <- function(rds_filename, learner, task) {
  m <- readRDS(file.path(model_root, learner, task, rds_filename))
  vals <- m$param_set$values
  prefix_re <- paste0("^", learner, "\\.")
  prefixed <- grep(prefix_re, names(vals), value = TRUE)
  if (length(prefixed) > 0) {
    out <- vals[prefixed]
    names(out) <- sub(prefix_re, "", names(out))
  } else {
    out <- vals
  }
  out[setdiff(names(out), reserved_keys)]
}

#### Build / load HP cache ####################################################

if (file.exists(hp_cache_location)) {
  message("Loading cached HPs from ", hp_cache_location)
  hp_cache <- readRDS(hp_cache_location)
} else {
  hp_cache <- list()
  for (task_name in names(preds_by_task)) {
    td <- preds_by_task[[task_name]]
    if (is.null(td$rds)) {
      stop("preds_by_task entry for ", task_name, " has no `rds` column. ",
           "Delete data/preds_cache.rds and rerun main.R to regenerate it.")
    }
    hp_cache[[task_name]] <- list()
    for (learner in unique(td$learner)) {
      in_family <- td$learner == learner
      n_fam <- sum(in_family)
      message(sprintf("  Loading %d %s models for %s ...", n_fam, learner, task_name))
      rdss <- td$rds[in_family]
      hp_list <- lapply(rdss, extract_hps_from_model, learner = learner, task = task_name)

      # Pad to common key set (some HPs may only exist for some models)
      all_keys <- unique(unlist(lapply(hp_list, names)))
      hp_df <- as.data.table(do.call(rbind, lapply(hp_list, function(hp) {
        row <- setNames(vector("list", length(all_keys)), all_keys)
        for (k in all_keys) row[[k]] <- if (k %in% names(hp)) hp[[k]] else NA
        row
      })))
      hp_cache[[task_name]][[learner]] <- hp_df
    }
  }
  saveRDS(hp_cache, hp_cache_location)
  message("Saved HP cache to ", hp_cache_location)
}

#### Plots ####################################################################

for (task_name in names(preds_by_task)) {
  td <- preds_by_task[[task_name]]
  for (learner in unique(td$learner)) {
    in_family <- td$learner == learner
    n_fam <- sum(in_family)
    if (n_fam < min_family_size) {
      message(task_name, " / ", learner, ": only ", n_fam,
              " models in RS -- skipping")
      next
    }

    fam_preds <- td$preds[in_family, , drop = FALSE]
    fam_mds1  <- cmdscale(dist(fam_preds), k = 1)[, 1]

    hp_df <- hp_cache[[task_name]][[learner]]
    if (is.null(hp_df) || nrow(hp_df) != n_fam) {
      message(task_name, " / ", learner,
              ": HP cache mismatch (have ", NROW(hp_df), " rows, expected ",
              n_fam, ") -- skipping")
      next
    }

    hp_cols <- names(hp_df)[vapply(names(hp_df), function(cn) {
      vals <- unlist(hp_df[[cn]])
      length(unique(vals[!is.na(vals)])) > 1
    }, logical(1))]

    if (length(hp_cols) == 0) {
      message(task_name, " / ", learner,
              ": no varying HPs in this family -- skipping")
      next
    }

    message(task_name, " / ", learner, ": plotting ", length(hp_cols),
            " HP(s) for ", n_fam, " models")

    n_per_row <- min(3, length(hp_cols))
    n_panel_rows <- ceiling(length(hp_cols) / n_per_row)
    pdf(file.path(hp_dir, sprintf("hp_vs_mds1_%s_%s.pdf", task_name, learner)),
        width = 4 * n_per_row, height = 4 * n_panel_rows)
    par(mfrow = c(n_panel_rows, n_per_row), mar = c(4.2, 4.2, 2.5, 1))
    for (hp in hp_cols) {
      val <- unlist(hp_df[[hp]])
      if (is.factor(val) || is.character(val) || is.logical(val)) {
        plot(factor(val), fam_mds1,
             xlab = hp, ylab = "MDS 1 (within-family)",
             main = sprintf("%s / %s", task_name, learner))
      } else {
        plot(val, fam_mds1,
             xlab = hp, ylab = "MDS 1 (within-family)",
             main = sprintf("%s / %s", task_name, learner),
             pch = 19, cex = 0.8)
      }
    }
    dev.off()
  }
}

message("Wrote HP trajectory plots to ", hp_dir, "/")
