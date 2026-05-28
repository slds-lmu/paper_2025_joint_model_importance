source("init.R")

#### settings #################################################################

run_models_merged_location <- "data/run_models_merged.rds"
preds_cache_location       <- "data/preds_cache.rds"
figures_dir                <- "figures"

# Rashomon set epsilon
RS_epsilon       <- 0.01
distance_metrics <- c("euclidean", "manhattan")

dir.create(figures_dir, showWarnings = FALSE)

#### prerequisites ############################################################

# performance and SVM kernel specification
perfor_and_kernel <- get_performance_and_SVMkernel(run_models_merged_location)
performance <- perfor_and_kernel$performance
kernel      <- perfor_and_kernel$kernel
rm(perfor_and_kernel)

# model indices of Rashomon set per task
RS <- get_RS(RS_epsilon, performance, vic)

#### compute predictions (batchtools, cached to disk) #########################
# Computing predictions for every Rashomon-set model is the expensive step.
# Once cached, downstream iteration on distances / clustering / plots needs
# no cluster access. Delete preds_cache_location to force a recompute.

if (file.exists(preds_cache_location)) {
  message("Loading cached predictions from ", preds_cache_location)
  preds_by_task <- readRDS(preds_cache_location)
} else {
  regr <- getRegistry(
    regpath      = "/media/external/ewaldf/paper_2025_joint_model_importance",
    make.default = TRUE
  )
  # writeable = TRUE (default) only in one window !!!!
  regr$cluster.functions <- makeClusterFunctionsSocket(ncpus = 30)

  addProblem("fromlist", fun = function(data, job, taskname) {
    task <- list.tasks[[taskname]]
    if (taskname == "bs") task <- fix_bs_task_for_featureimp(task)
    generateCanonicalDataSplits(task, ratio = 2 / 3, seed = 1)$validation
  })

  addAlgorithm("calculate_preds", fun = function(data, instance, job, learnername, model.no) {
    taskname <- job$pars$prob.pars$taskname
    model_path <- sprintf(
      "/media/external/rashomon/datafiles/%s/%s/samplemodel_%s_%s_%04d.rds",
      taskname, learnername, learnername, taskname, model.no
    )
    model <- readRDS(model_path)

    if (taskname == "bs") model <- fix_bs_model_for_predict(model, instance)

    if (instance$task_type == "classif") {
      if (model$predict_type != "prob") model$predict_type <- "prob"
      # binary tasks (gc, cs): one prob column suffices (the other is 1 - x)
      model$predict(instance)$prob[, 2]
    } else {
      model$predict(instance)$response
    }
  })

  pred_design <- rbindlist(
    lapply(names(RS), function(i) design[rn == i][RS[[i]]])
  )

  addExperiments(
    prob.designs = list(fromlist = data.table(taskname = pred_design$rn)),
    algo.designs = list(calculate_preds = pred_design[, -"rn"]),
    repls        = 1,
    combine      = "bind"
  )

  submitJobs(findErrors())
  submitJobs()
  waitForJobs()

  res <- unwrap(
    ijoin(
      getJobPars(),
      reduceResultsDataTable(fun = function(x) list(res = x))
    ),
    sep = "."
  )

  # one entry per task: matrix (row = model, col = validation observation)
  # plus the learner / model.no vectors so we can colour plots by learner.
  preds_by_task <- lapply(split(res, by = "taskname"), function(rres) {
    preds_mat <- do.call(rbind, rres$result.res)
    rownames(preds_mat) <- sprintf("%s_m%d", rres$learnername, rres$model.no)
    list(
      preds    = preds_mat,
      learner  = rres$learnername,
      model_no = rres$model.no
    )
  })

  saveRDS(preds_by_task, preds_cache_location)
  message("Saved predictions to ", preds_cache_location)
}

#### distances, clustering & visualisation ####################################
# For each task × distance metric, compute the pairwise behavioural distance
# between Rashomon-set models (rows = models, cols = validation observations)
# and visualise it three ways: hierarchical clustering, classical MDS,
# partitioning around medoids.

for (taskname in names(preds_by_task)) {
  td <- preds_by_task[[taskname]]
  learner_factor <- factor(td$learner)

  for (metric in distance_metrics) {
    dist_mat <- dist(td$preds, method = metric)

    # hierarchical clustering dendrogram
    pdf(file.path(figures_dir, sprintf("dendrogram_%s_%s.pdf", taskname, metric)),
        width = 9, height = 6)
    plot(hclust(dist_mat), hang = -1,
         main = sprintf("Hierarchical clustering — %s (%s)", taskname, metric),
         xlab = "model", sub = "", cex = 0.6)
    dev.off()

    # classical (metric) MDS, coloured by learner family
    mds <- cmdscale(dist_mat, k = 2, eig = TRUE)
    pdf(file.path(figures_dir, sprintf("mds_%s_%s.pdf", taskname, metric)),
        width = 7, height = 6)
    plot(mds$points, type = "n",
         xlab = "MDS 1", ylab = "MDS 2", asp = 1,
         main = sprintf("Classical MDS — %s (%s)", taskname, metric))
    points(mds$points,
           col = as.integer(learner_factor), pch = 19, cex = 0.9)
    legend("topright",
           legend = levels(learner_factor),
           col    = seq_len(nlevels(learner_factor)),
           pch    = 19, bty = "n")
    dev.off()

    # PAM, k = 2 placeholder — tune via silhouette in a later iteration
    pam_res <- cluster::pam(dist_mat, k = 2, metric = metric)
    pdf(file.path(figures_dir, sprintf("pam_%s_%s.pdf", taskname, metric)),
        width = 9, height = 5)
    plot(pam_res, main = sprintf("PAM clustering — %s (%s)", taskname, metric))
    dev.off()
  }
}

message("Wrote plots to ", figures_dir, "/")
