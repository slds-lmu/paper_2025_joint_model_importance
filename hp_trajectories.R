# Hyperparameter trajectories along the within-family MDS axis.
#
# Reads data/preds_cache.rds (written by main.R) and, for every
# (task, learner-family) with enough Rashomon-set members, projects each
# model onto a single 1D coordinate via within-family classical MDS, then
# plots that coordinate against each hyperparameter that varies across
# the family. One multi-panel PDF per (task, learner) into
# figures/hp_trajectories/.
#
# Assumption: in run_models$torun.samples[[learner]] filtered to a given
# taskname, the i-th row corresponds to model.no = i. The script prints
# a message() if this assumption looks violated for a family.

source("init.R")

preds_cache_location       <- "data/preds_cache.rds"
run_models_merged_location <- "data/run_models_merged.rds"
figures_dir                <- "figures"
hp_dir                     <- file.path(figures_dir, "hp_trajectories")
min_family_size            <- 5

if (!file.exists(preds_cache_location)) {
  stop("No predictions cache at ", preds_cache_location,
       " — run main.R first.")
}
preds_by_task <- readRDS(preds_cache_location)
run_models    <- readRDS(run_models_merged_location)

dir.create(hp_dir, showWarnings = FALSE, recursive = TRUE)

# Columns that are not hyperparameters in run_models$torun.samples[[learner]]
reserved_cols <- c("taskname", "score", "scores", "model.no", "id",
                   "iteration", "fold", "learnername", "row", "rn",
                   "batch.id", "trial.no", "repl")

# Loop variable name avoids the "taskname" column on hp_table — otherwise
# data.table's i-expression cannot disambiguate column vs. variable.
for (task_name in names(preds_by_task)) {
  td <- preds_by_task[[task_name]]
  for (learner in unique(td$learner)) {
    in_family <- td$learner == learner
    n_fam <- sum(in_family)
    if (n_fam < min_family_size) {
      message(task_name, " / ", learner, ": only ", n_fam,
              " models in RS — skipping")
      next
    }

    fam_preds    <- td$preds[in_family, , drop = FALSE]
    fam_mds1     <- cmdscale(dist(fam_preds), k = 1)[, 1]
    fam_model_no <- td$model_no[in_family]

    if (is.null(run_models$torun.samples[[learner]])) {
      message(task_name, " / ", learner,
              ": no torun.samples entry — skipping")
      next
    }
    hp_table <- as.data.table(run_models$torun.samples[[learner]])
    if (!"taskname" %in% colnames(hp_table)) {
      message(task_name, " / ", learner,
              ": no taskname column in torun.samples — skipping")
      next
    }
    hp_task <- hp_table[taskname == task_name]

    if (nrow(hp_task) < max(fam_model_no)) {
      message(task_name, " / ", learner, ": HP table has only ",
              nrow(hp_task), " rows but max model.no is ",
              max(fam_model_no), " — assumption broken, skipping")
      next
    }

    # The crucial assumption: model.no is a 1-based row index into hp_task
    hp_rows <- hp_task[fam_model_no]

    candidate_cols <- setdiff(colnames(hp_rows), reserved_cols)
    hp_cols <- candidate_cols[vapply(candidate_cols, function(cn) {
      length(unique(hp_rows[[cn]])) > 1
    }, logical(1))]

    if (length(hp_cols) == 0) {
      message(task_name, " / ", learner,
              ": no varying HP columns in this family — skipping")
      next
    }

    message(task_name, " / ", learner, ": plotting ", length(hp_cols),
            " HP(s) for ", n_fam, " models")

    n_per_row <- min(3, length(hp_cols))
    n_panel_rows <- ceiling(length(hp_cols) / n_per_row)
    pdf(file.path(hp_dir, sprintf("hp_vs_mds1_%s_%s.pdf", task_name, learner)),
        width = 4 * n_per_row, height = 4 * n_panel_rows)
    par(mfrow = c(n_panel_rows, n_per_row),
        mar = c(4.2, 4.2, 2.5, 1))
    for (hp in hp_cols) {
      val <- hp_rows[[hp]]
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
