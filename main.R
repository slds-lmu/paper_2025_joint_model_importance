source("init.R")

#### settings #################################################################

preds_cache_location <- "data/preds_cache.rds"
figures_dir          <- "figures"

# Rashomon set epsilon. The floor is an absolute fallback for tasks where
# best_score ~ 0 (e.g. cr, mk) -- without it best * epsilon collapses to
# zero and no model qualifies. For tasks with non-trivial best scores the
# relative term dominates and the floor is inactive.
RS_epsilon       <- 0.01
RS_epsilon_floor <- 0.005
distance_metrics <- c("euclidean", "manhattan")

dir.create(figures_dir, showWarnings = FALSE)

#### Rashomon set #############################################################
# Performance scores per (task, learner, model.no) come pre-computed in
# data/results_modelperformances.RData (loaded as `res_dt` by init.R).
# Build the nested list format get_RS() expects, then ask which models lie
# within epsilon of each task's best.

performance <- build_performance_from_res_dt(res_dt,
                                              tasks    = names(vic),
                                              learners = learner.keys)
RS <- get_RS(RS_epsilon, performance, vic, epsilon_floor = RS_epsilon_floor)

# Rashomon-set composition per task, useful sanity check
message("Learner set (auto-detected from VIC): ",
        paste(learner.keys, collapse = ", "))
message("Rashomon-set sizes (n_models per task):")
for (t in names(RS)) {
  message(sprintf("  %s: %d", t, length(RS[[t]])))
}

#### prediction extraction ####################################################
# Predictions on the held-out validation split are pre-computed by the
# upstream pred_mult.R pipeline in paper_2024_rashomon_set and live in
# data/results_preds_all_but_TreeFARMS.RData as a deeply-nested list
# preds[[task]][[learner]][[model.no]] of mlr3 Prediction objects.
# Here we slim that down to one numeric vector per RS member and cache
# the result so downstream iteration on plots does not have to reload
# the 425 MB predictions object.

extract_pred_vec <- function(pred_obj) {
  if (inherits(pred_obj, "PredictionClassif")) {
    # binary tasks (gc, cs): one probability column suffices (other = 1 - x)
    pred_obj$prob[, 2]
  } else {
    pred_obj$response
  }
}

if (file.exists(preds_cache_location)) {
  message("Loading cached predictions from ", preds_cache_location)
  preds_by_task <- readRDS(preds_cache_location)
} else {
  message("Loading predictions from data/results_preds_all_but_TreeFARMS.RData",
          " (this is ~425 MB and may take a moment)")
  load("data/results_preds_all_but_TreeFARMS.RData")  # provides `preds`

  preds_by_task <- lapply(names(RS), function(task_name) {
    rs_rows <- design[rn == task_name][RS[[task_name]]]
    if (nrow(rs_rows) == 0) return(NULL)
    pred_mat <- do.call(rbind, lapply(seq_len(nrow(rs_rows)), function(i) {
      extract_pred_vec(preds[[task_name]][[rs_rows$learnername[i]]][[rs_rows$model.no[i]]])
    }))
    rownames(pred_mat) <- sprintf("%s_m%d", rs_rows$learnername, rs_rows$model.no)
    list(
      preds    = pred_mat,
      learner  = rs_rows$learnername,
      model_no = rs_rows$model.no,
      rds      = rs_rows$rds
    )
  })
  names(preds_by_task) <- names(RS)
  preds_by_task <- preds_by_task[!sapply(preds_by_task, is.null)]

  saveRDS(preds_by_task, preds_cache_location)
  message("Saved predictions to ", preds_cache_location)
  rm(preds)  # free the 425 MB
  invisible(gc())
}

#### distances, clustering & visualisation ####################################
# For each task x distance metric, compute the pairwise behavioural distance
# between Rashomon-set models (rows = models, cols = validation observations)
# and visualise it three ways: hierarchical clustering, classical MDS,
# partitioning around medoids.

for (task_name in names(preds_by_task)) {
  td <- preds_by_task[[task_name]]
  if (nrow(td$preds) < 3) {
    message(task_name, ": only ", nrow(td$preds),
            " RS member(s) -- need >= 3 for MDS(k=2) / PAM(k=2), skipping")
    next
  }
  learner_factor <- factor(td$learner)

  for (metric in distance_metrics) {
    dist_mat <- dist(td$preds, method = metric)

    # hierarchical clustering dendrogram
    pdf(file.path(figures_dir, sprintf("dendrogram_%s_%s.pdf", task_name, metric)),
        width = 9, height = 6)
    plot(hclust(dist_mat), hang = -1,
         main = sprintf("Hierarchical clustering -- %s (%s)", task_name, metric),
         xlab = "model", sub = "", cex = 0.6)
    dev.off()

    # classical (metric) MDS, coloured by learner family
    mds <- cmdscale(dist_mat, k = 2, eig = TRUE)
    pdf(file.path(figures_dir, sprintf("mds_%s_%s.pdf", task_name, metric)),
        width = 7, height = 6)
    plot(mds$points, type = "n",
         xlab = "MDS 1", ylab = "MDS 2", asp = 1,
         main = sprintf("Classical MDS -- %s (%s)", task_name, metric))
    points(mds$points,
           col = as.integer(learner_factor), pch = 19, cex = 0.9)
    legend("topright",
           legend = levels(learner_factor),
           col    = seq_len(nlevels(learner_factor)),
           pch    = 19, bty = "n")
    dev.off()

    # PAM, k = 2 placeholder -- tune via silhouette in a later iteration
    pam_res <- cluster::pam(dist_mat, k = 2, metric = metric)
    pdf(file.path(figures_dir, sprintf("pam_%s_%s.pdf", task_name, metric)),
        width = 9, height = 5)
    plot(pam_res, main = sprintf("PAM clustering -- %s (%s)", task_name, metric))
    dev.off()
  }
}

message("Wrote plots to ", figures_dir, "/")
